import os
import pandas as pd
from Bio import SeqIO
import multiprocessing
from panCG.lib.base import BasePan
from panCG.lib.pangeneIndex import GeneAnchor
from panCG.lib.base import ParallelScheduler


def run_geneMapMerge(logger, config, workDir, Reference, args_dict):
    logger.info("---------------------------------- step 1. Start mergeGeneMap ... -----------------------------------")
    logger.info("args: {}".format(args_dict))
    configData = BasePan.read_yaml(config)

    species_tree = configData["species_tree"]
    speciesList = BasePan.get_species_li(species_tree, Reference)
    logger.info("speciesList: {}".format(speciesList))
    logger.info("species num: {}".format(len(speciesList)))

    # Check the configuration
    logger.info("Check the configuration ...")
    all_bed_pep_name_li = []
    for spe in speciesList:
        bed_file = configData["species"][spe]["longest_pep_bed"]
        fa_file = configData["species"][spe]["longest_pep_fasta"]
        fa_pep_name_li = []
        for record in SeqIO.parse(fa_file, "fasta"):
            fa_pep_name_li.append(record.id)
        bed_df = pd.read_csv(bed_file, sep="\t", header=None)
        bed_pep_name_li = bed_df[3].tolist()
        if not all(i in bed_pep_name_li for i in fa_pep_name_li):
            raise Exception(f"The bed and file ids in {spe} are inconsistent")
        all_bed_pep_name_li.extend(bed_pep_name_li)
    if len(set(all_bed_pep_name_li)) != len(all_bed_pep_name_li):
        raise Exception("The protein sequence ID of the input species is duplicated")
    logger.info("Configuration file is available")

    # Add class attributes to GeneAnchor class
    GeneAnchor.set_diamond_run_threads(args_dict["diamond_threads"])
    GeneAnchor.set_jcvi_run_threads(args_dict["jcvi_threads"])
    GeneAnchor.set_identity_threshold(args_dict["identity"])
    GeneAnchor.set_E_value(args_dict["evalue"])

    diamondDir = os.path.join(workDir, "diamondDir")
    JCVIDir = os.path.join(workDir, "JCVIDir")
    MergeDir = os.path.join(workDir, "Merge_JCVI_diamond")
    BasePan.pymkdir(workDir)
    BasePan.pymkdir(diamondDir)
    BasePan.pymkdir(JCVIDir)
    BasePan.pymkdir(MergeDir)

    # Create instances in batches
    GeneAnchorList = []
    for i in range(len(speciesList)):
        for j in range(i, len(speciesList)):
            que, ref = speciesList[i], speciesList[j]
            queProt = configData["species"][que]["longest_pep_fasta"]
            refProt = configData["species"][ref]["longest_pep_fasta"]
            queLongestProtBed = configData["species"][que]["longest_pep_bed"]
            refLongestProtBed = configData["species"][ref]["longest_pep_bed"]
            geneAnchor = GeneAnchor(que, ref, queProt, refProt, queLongestProtBed, refLongestProtBed, workDir)
            GeneAnchorList.append(geneAnchor)

    # 2. parallel diamond
    # diamond align
    diamond_ok_file = os.path.join(diamondDir, "diamond.ok")
    if not os.path.exists(diamond_ok_file):
        logger.info("parallel diamond start ...")
        diamond_threads = args_dict["diamond_parallel"]
        total_tasks = len(GeneAnchorList)
        logger.info("total_tasks: {}".format(total_tasks))
        logger.info("threads: {}".format(diamond_threads))
        completed_tasks = multiprocessing.Manager().Value('i', 0)
        results = []
        pool = multiprocessing.Pool(processes=diamond_threads)
        progress_block = diamond_threads
        parallelScheduler = ParallelScheduler(logger, total_tasks, completed_tasks, progress_block)
        for gAnchor in GeneAnchorList:
            outDir = os.path.join(diamondDir, f"{gAnchor.que}.{gAnchor.ref}")
            BasePan.pymkdir(outDir)
            track_progress = parallelScheduler.make_call_back()
            error_callback = parallelScheduler.make_error_callback()
            result = pool.apply_async(GeneAnchor.multi_RunDiamond, args=(gAnchor, outDir,),
                                      callback=track_progress, error_callback=error_callback)
            results.append(result)
        pool.close()
        pool.join()
        del total_tasks, completed_tasks, progress_block, parallelScheduler, results
        logger.info("diamond align Finish !!!")
        BasePan.touch_file(diamond_ok_file)
    else:
        logger.info("jump diamond ...")
    del diamond_ok_file

    # diamond filter by identity
    diamond_filter_ok_file = os.path.join(diamondDir, "diamond_filter.ok")
    if not os.path.exists(diamond_filter_ok_file):
        logger.info("diamond filter start ...")
        diamond_threads = args_dict["diamond_parallel"]
        total_tasks = len(GeneAnchorList)
        logger.info("total_tasks: {}".format(total_tasks))
        logger.info("threads: {}".format(diamond_threads))
        completed_tasks = multiprocessing.Manager().Value('i', 0)
        results = []
        pool = multiprocessing.Pool(processes=diamond_threads)
        progress_block = diamond_threads
        parallelScheduler = ParallelScheduler(logger, total_tasks, completed_tasks, progress_block)
        for gAnchor in GeneAnchorList:
            outDir = os.path.join(diamondDir, f"{gAnchor.que}.{gAnchor.ref}")
            BasePan.pymkdir(outDir)
            track_progress = parallelScheduler.make_call_back()
            error_callback = parallelScheduler.make_error_callback()
            result = pool.apply_async(GeneAnchor.multi_FilterDiamond, args=(gAnchor, outDir,),
                                      callback=track_progress, error_callback=error_callback)
            results.append(result)
        pool.close()
        pool.join()
        del total_tasks, completed_tasks, progress_block, parallelScheduler, results
        logger.info("diamond filter Finish ...")
        BasePan.touch_file(diamond_filter_ok_file)
    else:
        logger.info("jump diamond filter ...")
    del diamond_filter_ok_file

    # 2. parallel JCVI
    # JCVI align
    jcvi_ok_file = os.path.join(JCVIDir, "jcvi.ok")
    if not os.path.exists(jcvi_ok_file):
        logger.info("JCVI align start ...")
        jcvi_threads = args_dict["jcvi_parallel"]
        total_tasks = len(GeneAnchorList)
        logger.info("total_tasks: {}".format(total_tasks))
        logger.info("threads: {}".format(jcvi_threads))
        completed_tasks = multiprocessing.Manager().Value('i', 0)
        results = []
        pool = multiprocessing.Pool(processes=jcvi_threads)
        progress_block = jcvi_threads
        parallelScheduler = ParallelScheduler(logger, total_tasks, completed_tasks, progress_block)
        for gAnchor in GeneAnchorList:
            outDir = os.path.join(JCVIDir, f"{gAnchor.que}.{gAnchor.ref}")
            BasePan.pymkdir(outDir)
            track_progress = parallelScheduler.make_call_back()
            error_callback = parallelScheduler.make_error_callback()
            result = pool.apply_async(GeneAnchor.multi_RunJCVI, args=(gAnchor, outDir,),
                                      callback=track_progress, error_callback=error_callback)
            results.append(result)
        pool.close()
        pool.join()
        del total_tasks, completed_tasks, progress_block, parallelScheduler, results
        logger.info("JCVI align Finish")
        BasePan.touch_file(jcvi_ok_file)
    else:
        logger.info("jump JCVI align ...")
    del jcvi_ok_file

    # 3. get Two-way anchor
    logger.info("get Two-way anchor start ...")
    TwoWay_threads = args_dict["diamond_parallel"]
    total_tasks = len(GeneAnchorList)
    logger.info("total_tasks: {}".format(total_tasks))
    logger.info("threads: {}".format(TwoWay_threads))
    completed_tasks = multiprocessing.Manager().Value('i', 0)
    results = []
    pool = multiprocessing.Pool(processes=TwoWay_threads)
    progress_block = TwoWay_threads
    parallelScheduler = ParallelScheduler(logger, total_tasks, completed_tasks, progress_block)
    for gAnchor in GeneAnchorList:
        if gAnchor.que == gAnchor.ref:
            continue
        track_progress = parallelScheduler.make_call_back()
        error_callback = parallelScheduler.make_error_callback()
        result = pool.apply_async(GeneAnchor.multi_get_TowWey, args=(gAnchor,),
                                  callback=track_progress, error_callback=error_callback)
        results.append(result)
    pool.close()
    pool.join()
    del total_tasks, completed_tasks, progress_block, parallelScheduler, results
    logger.info("get Two-way anchor Finish")

    Species_pair_li = []
    for que in speciesList:
        for ref in speciesList:
            Species_pair_li.append((que, ref))
            # if que != ref:
            #     Species_pair_li.append((que, ref))

    # get anchor
    merge_tasks = []
    for i in Species_pair_li:
        que, ref = i[0], i[1]
        # diamond
        diamond_blastp_file = os.path.join(diamondDir, f"{que}.{ref}", f"{que}.{ref}.diamond_blastp.filter.fmt6.txt")
        diamond_blastp_file_anchor = os.path.join(diamondDir, f"{que}.{ref}", f"{que}.{ref}.diamond_blastp.filter.anchors")
        GeneAnchor.get_anchor(diamond_blastp_file, diamond_blastp_file_anchor)
        # jcvi
        anchors_file = os.path.join(JCVIDir, f"{que}.{ref}", f"{que}.{ref}.anchors")
        anchors_file_anchor = os.path.join(JCVIDir, f"{que}.{ref}", f"{que}.{ref}.filter.anchors")
        GeneAnchor.get_anchor(anchors_file, anchors_file_anchor)
        lifted_anchors_file = os.path.join(JCVIDir, f"{que}.{ref}", f"{que}.{ref}.lifted.anchors")
        lifted_anchors_file_anchor = os.path.join(JCVIDir, f"{que}.{ref}", f"{que}.{ref}.filter.lifted.anchors")
        GeneAnchor.get_anchor(lifted_anchors_file, lifted_anchors_file_anchor)

        BasePan.pymkdir(os.path.join(MergeDir, f"{que}.{ref}"))
        merge_out_file = os.path.join(MergeDir, f"{que}.{ref}", f"{que}.{ref}.merge.anchors")
        merge_tasks.append((merge_out_file, diamond_blastp_file_anchor, lifted_anchors_file_anchor))

    # merge
    logger.info("merge JCVI and diamond ...")
    merge_threads = args_dict["diamond_parallel"]
    total_tasks = len(Species_pair_li)
    logger.info("total_tasks: {}".format(total_tasks))
    logger.info("threads: {}".format(merge_threads))
    completed_tasks = multiprocessing.Manager().Value('i', 0)
    results = []
    pool = multiprocessing.Pool(processes=merge_threads)
    progress_block = merge_threads
    parallelScheduler = ParallelScheduler(logger, total_tasks, completed_tasks, progress_block)
    for i in merge_tasks:
        track_progress = parallelScheduler.make_call_back()
        error_callback = parallelScheduler.make_error_callback()
        result = pool.apply_async(BasePan.merge_files, args=i,
                                  callback=track_progress, error_callback=error_callback)
        results.append(result)
    pool.close()
    pool.join()
    del total_tasks, completed_tasks, progress_block, parallelScheduler, results
    logger.info("merge Finish")


