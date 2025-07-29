import os
import pandas as pd
import multiprocessing
import pyBigWig
import numpy as np
import pyranges as pr

from panCG.lib.base import ParallelScheduler
from panCG.lib.base import BasePan
from panCG.lib.base import split_base_group_data
from panCG.lib.cns_recall import HalLiftoverMapCNSs
from panCG.lib.cns_recall import ReCnsIndexMerge
from panCG.lib.cns_recall import CeCds
from panCG.lib.cns_recall import TripleCnsIndexMerge


def getCDS_from_GffFile(gffFile, cdsFile):
    gff_df = pd.read_csv(gffFile, sep="\t", header=None, comment='#', skip_blank_lines=True)
    gff_df = gff_df[gff_df[2] == "CDS"]
    gff_df[3] = gff_df[3] - 1
    gff_df = gff_df[[0, 3, 4]]
    gff_df.columns = ["Chromosome", "Start", "End"]
    gff_pr = pr.PyRanges(gff_df)
    gff_pr = gff_pr.sort().merge()
    gff_data = gff_pr.df
    gff_data.to_csv(cdsFile, sep="\t", header=False, index=False)


def run_Recall(data, speciesList, workDir, all_halLiftover_dict):
    datas_li = []
    grouped = data.groupby('Group', sort=True)
    for group_name, group_data in grouped:
        HalLiftoverMapCNSser = HalLiftoverMapCNSs(group_data, all_halLiftover_dict, speciesList, workDir)
        datas_li.append(HalLiftoverMapCNSser.run())
    return pd.concat(datas_li, ignore_index=True, axis=0)


def add_bw_score(merge_bed_file, bw_file, out_file):
    """
    This function is used to add two columns to the bed file, namely averageBwScore and effective_len
    Args:
        merge_bed_file: A bed file with four columns, the fourth column is the name of the cns
        bw_file: phastCons bw file
        out_file: A bed-like file with 6 columns
    Returns:
        None
    """
    bw = pyBigWig.open(bw_file, "r")
    txt = ""
    with open(merge_bed_file, "r") as f:
        for line in f:
            line_li = line.strip().split("\t")
            chrID, start, end = line_li[0], int(line_li[1]), int(line_li[2])
            averageBwScore_tmp = bw.stats(chrID, start, end)[0]
            averageBwScore = averageBwScore_tmp if averageBwScore_tmp is not None else 0  # For None
            effective_len = np.count_nonzero(~np.isnan(bw.values(chrID, start, end)))
            line_li.extend([str(averageBwScore), str(effective_len)])
            txt += "\t".join(line_li) + "\n"
    bw.close()
    fo = open(out_file, "w")
    fo.write(txt)
    fo.close()


def parse_merge_bed(merge_bed_file):
    out_dict = {}
    with open(merge_bed_file, "r") as f:
        for line in f:
            line_li = line.strip().split("\t")
            map_region = [(line_li[0], int(line_li[1]), int(line_li[2])), (float(line_li[4]), int(line_li[5]))]
            out_dict.setdefault(line_li[3], []).append(map_region)
    return out_dict


def multi_ReCnsIndexMerge(data, speciesList):
    Grouped = data.groupby('Group', sort=False)
    group_data_li = []
    for _, group_data in Grouped:
        if group_data.shape[0] == 1:
            group_data_li.append(group_data)
        else:
            ReCnsIndexMerger = ReCnsIndexMerge(group_data, speciesList)
            ReCnsIndexMerger.indexMerge()
            group_data_li.append(ReCnsIndexMerger.group_data)
    return pd.concat(group_data_li, axis=0)


class ReplaceCdsCe:
    @staticmethod
    def replace_cds_ce(value, set_li):
        if value == ".":
            return value
        else:
            out_li = []
            values_li = value.split(",")
            for i in values_li:
                if "(CEs)" in i:
                    if i.replace("(CEs)", "") in set_li:
                        out_li.append(i.replace("(CEs)", "(recall_cds)"))
                    else:
                        out_li.append(i.replace("(CEs)", "(recall_cns)"))
                else:
                    out_li.append(i)
            return ",".join(out_li)


def multi_TripleCnsIndexMerge(data, speciesList):
    Grouped = data.groupby('Group', sort=False)
    group_data_li = []
    for _, group_data in Grouped:
        if group_data.shape[0] == 1:
            group_data_li.append(group_data)
        else:
            TripleCnsIndexMerger = TripleCnsIndexMerge(group_data, speciesList)
            TripleCnsIndexMerger.indexMerge()
            group_data_li.append(TripleCnsIndexMerger.group_data)
    return pd.concat(group_data_li, axis=0)


def run_cnsRecall(logger, config, workDir, Reference, args_dict):
    """

    Args:
        args_dict:
        logger: logger
        config: config
        workDir: workDir
        Reference: Reference

    Returns:

    """
    logger.info("---------------------------------- step 5. Start run_cnsRecall ... ----------------------------------")
    RefIndexDir = os.path.join(workDir, "Ref_{}_IndexDir".format(Reference))
    cnsIndexMergeFile = os.path.join(RefIndexDir, "Ref.{}.cnsIndexMerge.csv".format(Reference))
    cnsIndexMerge_data = pd.read_table(cnsIndexMergeFile)

    configData = BasePan.read_yaml(config)
    species_Dict = configData["species"]
    halLiftoverDir = os.path.join(workDir, "halLiftoverDir")

    species_tree = configData["species_tree"]
    speciesList = BasePan.get_species_li(species_tree, Reference)
    species_pairs = []
    for que in speciesList:
        for ref in speciesList:
            if que != ref:
                species_pairs.append((que, ref))

    # Calculate the conservative score of the combined HalLiftover results
    logger.info("Start pyBigWig ...")
    pyBigWig_threads = args_dict["recall_threads"]
    total_tasks = len(species_pairs)
    logger.info("total_tasks: {}".format(total_tasks))
    logger.info("threads: {}".format(pyBigWig_threads))
    completed_tasks = multiprocessing.Manager().Value('i', 0)
    progress_block = pyBigWig_threads
    parallelScheduler = ParallelScheduler(logger, total_tasks, completed_tasks, progress_block)
    pool = multiprocessing.Pool(processes=pyBigWig_threads)
    results = []
    for species_pair in species_pairs:
        que, ref = species_pair[0], species_pair[1]
        merge_file = os.path.join(halLiftoverDir, f"{que}.{ref}", f"{que}.{ref}.merge.bed")
        bw_file = configData["species"][ref]["bwFile"]
        merge_bw_bed_file = os.path.join(halLiftoverDir, f"{que}.{ref}", f"{que}.{ref}.merge.bw.bed")
        track_progress = parallelScheduler.make_call_back()
        error_callback = parallelScheduler.make_error_callback()
        result = pool.apply_async(add_bw_score,
                                  args=(merge_file, bw_file, merge_bw_bed_file,),
                                  callback=track_progress, error_callback=error_callback)
        results.append(result)
    pool.close()
    pool.join()
    logger.info("End pyBigWig.")

    # Store the result of the previous step into a dictionary
    logger.info("Start create all_halLiftover_dict ...")
    all_halLiftover_dict = {}
    for species_pair in species_pairs:
        que, ref = species_pair[0], species_pair[1]
        merge_bed_file = os.path.join(halLiftoverDir, f"{que}.{ref}", f"{que}.{ref}.merge.bw.bed")
        out_dict = parse_merge_bed(merge_bed_file)
        all_halLiftover_dict[f"{que}.{ref}"] = out_dict
    logger.info("End create all_halLiftover_dict ...")

    # Add class attribute to HalLiftoverMapCNSs class
    HalLiftoverMapCNSs.set_AverageBwScoreThreshold(args_dict["aver_bw_score_threshold"])
    HalLiftoverMapCNSs.set_effective_rate(args_dict["recall_effective_rate"])
    HalLiftoverMapCNSs.set_cns_rate(args_dict["recall_cns_rate"])

    # Generate cds bed file
    logger.info("Start Generate cds bed file ...")
    cdsDir = os.path.join(workDir, "cdsDir")
    BasePan.pymkdir(cdsDir)
    logger.info("speciesList: {}".format(speciesList))
    logger.info("len(speciesList): {}".format(len(speciesList)))
    for specie in speciesList:
        cdsFile = os.path.join(cdsDir, f"{specie}.sort.cds.bed")
        if not os.path.exists(cdsFile):
            gffFile = species_Dict[specie]["gffFile"]
            getCDS_from_GffFile(gffFile, cdsFile)
    logger.info("End Generate cds bed file ...")

    # run recall
    # This step requires a large amount of memory.
    # A single Brassica thread requires 15Gb, but the running time is short and there is no need to parallelize it.
    logger.info("Start recall ...")
    recall_file = os.path.join(RefIndexDir, "Ref.{}.recallCEs.csv".format(Reference))
    recall_data = run_Recall(cnsIndexMerge_data, speciesList, workDir, all_halLiftover_dict)
    recall_data.to_csv(recall_file, sep='\t', index=False)
    logger.info("End recall ...")

    # Index after merge recall
    logger.info("Start ReCnsIndexMerge ...")
    block_size = args_dict["merge_chunk_size"]
    blockData_li = split_base_group_data(recall_data, 'Group', block_size)
    total_tasks = len(blockData_li)
    threads = args_dict["merge_threads"]
    logger.info("total_tasks: {}".format(total_tasks))
    logger.info("threads: {}".format(threads))
    completed_tasks = multiprocessing.Manager().Value('i', 0)
    progress_block = threads
    parallelScheduler = ParallelScheduler(logger, total_tasks, completed_tasks, progress_block)
    pool = multiprocessing.Pool(processes=threads)
    results = []
    for blockData in blockData_li:
        track_progress = parallelScheduler.make_call_back()
        error_callback = parallelScheduler.make_error_callback()
        result = pool.apply_async(multi_ReCnsIndexMerge, args=(blockData, speciesList,),
                                  callback=track_progress, error_callback=error_callback)
        results.append(result)
    pool.close()
    pool.join()
    ReCnsIndexMergeData = pd.concat([i.get() for i in results], axis=0)
    ReCnsIndexMergeData['Index'] = ReCnsIndexMergeData.groupby('Group').cumcount() + 1
    ReCnsIndexMergeData['Index'] = ReCnsIndexMergeData['Group'] + '.' + ReCnsIndexMergeData['Index'].astype(str)
    ReCnsIndexMerge_file = os.path.join(RefIndexDir, "Ref.{}.ReCnsIndexMerge.csv".format(Reference))
    ReCnsIndexMergeData.to_csv(ReCnsIndexMerge_file, sep='\t', index=False)
    logger.info("End ReCnsIndexMerge ...")

    # Get CE that overlaps with cds
    CEsDir = os.path.join(workDir, "CEsDir")
    BasePan.pymkdir(CEsDir)
    for species in speciesList:
        cds_file = os.path.join(workDir, "cdsDir", f"{species}.sort.cds.bed")
        out_file = os.path.join(CEsDir, f"{species}.recall_cds.bed")
        CeCdser = CeCds(ReCnsIndexMergeData, species, cds_file, out_file)
        CeCdser.get_ces()

    # Replace CEs with recall_cns or recall_cds
    CEsDir = os.path.join(workDir, "CEsDir")
    for species in speciesList:
        cds_ce_file = os.path.join(CEsDir, f"{species}.recall_cds.bed")
        cds_ce_data = pd.read_table(cds_ce_file, header=None, names=["chrID", "start", "end", "cds_ce_id"])
        cds_ce_set = set(cds_ce_data["cds_ce_id"].tolist())
        ReCnsIndexMergeData[species] = (ReCnsIndexMergeData[species]
                                        .apply(lambda x: ReplaceCdsCe.replace_cds_ce(x, cds_ce_set)))

    # Merge indexes with a large number of overlaps in the same group
    logger.info("Start TripleCnsIndexMerge ...")
    block_size = args_dict["recall_chunk_size"]
    blockData_li = split_base_group_data(ReCnsIndexMergeData, 'Group', block_size)
    total_tasks = len(blockData_li)
    threads = args_dict["merge_threads"]
    logger.info("total_tasks: {}".format(total_tasks))
    logger.info("threads: {}".format(threads))
    completed_tasks = multiprocessing.Manager().Value('i', 0)
    progress_block = args_dict["merge_threads"]
    parallelScheduler = ParallelScheduler(logger, total_tasks, completed_tasks, progress_block)
    pool = multiprocessing.Pool(processes=threads)
    results = []
    for blockData in blockData_li:
        track_progress = parallelScheduler.make_call_back()
        error_callback = parallelScheduler.make_error_callback()
        result = pool.apply_async(multi_TripleCnsIndexMerge,
                                  args=(blockData, speciesList,),
                                  callback=track_progress,
                                  error_callback=error_callback)
        results.append(result)
    pool.close()
    pool.join()
    TripleCnsIndexMergeData = pd.concat([i.get() for i in results], axis=0)
    TripleCnsIndexMergeData['Index'] = TripleCnsIndexMergeData.groupby('Group').cumcount() + 1
    TripleCnsIndexMergeData['Index'] = TripleCnsIndexMergeData['Group'] + '.' + TripleCnsIndexMergeData['Index'].astype(str)
    TripleCnsIndexMerge_file = os.path.join(RefIndexDir, "Ref.{}.TripleCnsIndexMerge.csv".format(Reference))
    TripleCnsIndexMergeData.to_csv(TripleCnsIndexMerge_file, sep='\t', index=False)
    logger.info("End TripleCnsIndexMerge ...")

    # Merge CEs of the same type that overlap in the same cell
    logger.info("Merge CEs of the same type that overlap in the same cell ...")
    for species in speciesList:
        TripleCnsIndexMergeData[species] = TripleCnsIndexMergeData[species].apply(lambda x: TripleCnsIndexMerge.merge_overlap(x))

    cnsRecall_file = os.path.join(RefIndexDir, "Ref.{}.recall.csv".format(Reference))
    TripleCnsIndexMergeData.to_csv(cnsRecall_file, sep='\t', index=False)

    logger.info("End run_cnsRecall ...")

