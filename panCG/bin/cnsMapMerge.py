"""
This script first obtains the correspondence between the CNS in ref between
the CNS in que through the `halLiftover` command in Cactus,
and then uses this to filter the results of blastn and JCVI.
Then get the corresponding relationship between ref and que's final CNS.
"""

import os
import multiprocessing
from panCG.lib.base import BasePan
from panCG.lib.pancns import CNSsAnchor
from panCG.lib.base import ParallelScheduler
from panCG.lib.pancns import run_rescue_cns_anchor


def run_cnsMapMerge(logger, config, workDir, geneWorkDir, reference, args_dict):
    logger.info("----------------------------------- step 1. Start cnsMapMerge ... ------------------------------------")
    logger.info("args: {}".format(args_dict))
    configData = BasePan.read_yaml(config)

    # Add class attributes to CNSsAnchor class
    CNSsAnchor.modify_blastn_evalue(args_dict["blastn_evalue"])
    CNSsAnchor.modify_max_gap(args_dict["map_gap"])
    CNSsAnchor.modify_threads(args_dict["blastn_threads"], args_dict["jcvi_threads"])
    CNSsAnchor.modify_overlap_rate(args_dict["halLiftover_rate"])

    species_Dict = configData["species"]
    HalFile = configData["HalFile"]
    species_tree = configData["species_tree"]
    speciesList = BasePan.get_species_li(species_tree, reference)

    halLiftoverDir = os.path.join(workDir, "halLiftoverDir")
    blastnDir = os.path.join(workDir, "blastnDir")
    JCVIDir = os.path.join(workDir, "JCVIDir")
    BasePan.pymkdir(workDir)
    BasePan.pymkdir(halLiftoverDir)
    BasePan.pymkdir(blastnDir)
    BasePan.pymkdir(JCVIDir)

    # Create instance for halLiftover
    Pair_CNSsAnchorList = []
    for que in speciesList:
        for ref in speciesList:
            if que == ref:
                continue
            queCNSsBedFile = species_Dict[que]["CNS_bed"]
            refCNSsBedFile = species_Dict[ref]["CNS_bed"]
            queGenomeFile = species_Dict[que]["GenomeFile"]
            refGenomeFile = species_Dict[ref]["GenomeFile"]
            cnsAnchor = CNSsAnchor(HalFile, que, ref, queCNSsBedFile, refCNSsBedFile,
                                   queGenomeFile, refGenomeFile, workDir)
            Pair_CNSsAnchorList.append(cnsAnchor)

    # 1. parallel halLiftover
    logger.info("parallel halLiftover start ...")
    halLiftover_threads = args_dict["halLiftover_parallel"]
    total_tasks = len(Pair_CNSsAnchorList)
    logger.info("total_tasks: {}".format(total_tasks))
    logger.info("threads: {}".format(halLiftover_threads))
    completed_tasks = multiprocessing.Manager().Value('i', 0)
    results = []
    pool = multiprocessing.Pool(processes=halLiftover_threads)
    progress_block = halLiftover_threads
    parallelScheduler = ParallelScheduler(logger, total_tasks, completed_tasks, progress_block)
    for cnsAnchor in Pair_CNSsAnchorList:
        halLiftoverOutDir = os.path.join(halLiftoverDir, f"{cnsAnchor.que}.{cnsAnchor.ref}")
        BasePan.pymkdir(halLiftoverOutDir)
        track_progress = parallelScheduler.make_call_back()
        error_callback = parallelScheduler.make_error_callback()
        result = pool.apply_async(CNSsAnchor.multi_halLiftover, args=(cnsAnchor, halLiftoverOutDir,),
                                  callback=track_progress, error_callback=error_callback)
        results.append(result)
    pool.close()
    pool.join()
    del total_tasks, completed_tasks, progress_block, parallelScheduler, results
    logger.info("parallel halLiftover Finish")

    # Create instance for blastn
    CNSsAnchorList = []
    for i in range(len(speciesList)):
        for j in range(i + 1, len(speciesList)):
            que, ref = speciesList[i], speciesList[j]
            queCNSsBedFile = species_Dict[que]["CNS_bed"]
            refCNSsBedFile = species_Dict[ref]["CNS_bed"]
            queGenomeFile = species_Dict[que]["GenomeFile"]
            refGenomeFile = species_Dict[ref]["GenomeFile"]
            cnsAnchor = CNSsAnchor(HalFile, que, ref, queCNSsBedFile, refCNSsBedFile,
                                   queGenomeFile, refGenomeFile, workDir)
            CNSsAnchorList.append(cnsAnchor)

    # 2. parallel blastn
    logger.info("parallel blastn start ...")
    blastn_threads = args_dict["blastn_parallel"]
    total_tasks = len(CNSsAnchorList)
    logger.info("total_tasks: {}".format(total_tasks))
    logger.info("threads: {}".format(blastn_threads))
    completed_tasks = multiprocessing.Manager().Value('i', 0)
    results = []
    pool = multiprocessing.Pool(processes=blastn_threads)
    progress_block = blastn_threads
    parallelScheduler = ParallelScheduler(logger, total_tasks, completed_tasks, progress_block)
    for cnsAnchor in CNSsAnchorList:
        outDir = os.path.join(blastnDir, f"{cnsAnchor.que}.{cnsAnchor.ref}")
        BasePan.pymkdir(outDir)
        track_progress = parallelScheduler.make_call_back()
        error_callback = parallelScheduler.make_error_callback()
        result = pool.apply_async(CNSsAnchor.multi_RunBlastn, args=(cnsAnchor, outDir,),
                                  callback=track_progress, error_callback=error_callback)
        results.append(result)
    pool.close()
    pool.join()
    del total_tasks, completed_tasks, progress_block, parallelScheduler, results
    logger.info("parallel blastn Finish")

    # 3. parallel blastn halLiftover filter
    logger.info("parallel blastn halLiftover filter start ...")
    blastn_filter_threads = args_dict["blastn_parallel"]
    total_tasks = len(CNSsAnchorList)
    logger.info("total_tasks: {}".format(total_tasks))
    logger.info("threads: {}".format(blastn_filter_threads))
    completed_tasks = multiprocessing.Manager().Value('i', 0)
    results = []
    pool = multiprocessing.Pool(processes=blastn_filter_threads)
    progress_block = blastn_filter_threads
    parallelScheduler = ParallelScheduler(logger, total_tasks, completed_tasks, progress_block)
    for cnsAnchor in CNSsAnchorList:
        outDir = os.path.join(blastnDir, f"{cnsAnchor.que}.{cnsAnchor.ref}")
        track_progress = parallelScheduler.make_call_back()
        error_callback = parallelScheduler.make_error_callback()
        result = pool.apply_async(CNSsAnchor.multi_halLiftoverFilter_blastn, args=(cnsAnchor, outDir, workDir,),
                                  callback=track_progress, error_callback=error_callback)
        results.append(result)
    pool.close()
    pool.join()
    del total_tasks, completed_tasks, progress_block, parallelScheduler, results
    logger.info("parallel blastn halLiftover filter Finish")

    # 4. parallel jcvi align
    logger.info("parallel jcvi align start ...")
    JCVI_threads = args_dict["jcvi_parallel"]
    total_tasks = len(CNSsAnchorList)
    logger.info("total_tasks: {}".format(total_tasks))
    logger.info("threads: {}".format(JCVI_threads))
    completed_tasks = multiprocessing.Manager().Value('i', 0)
    results = []
    pool = multiprocessing.Pool(processes=JCVI_threads)
    progress_block = JCVI_threads
    parallelScheduler = ParallelScheduler(logger, total_tasks, completed_tasks, progress_block)
    for cnsAnchor in CNSsAnchorList:
        outDir = os.path.join(JCVIDir, f"{cnsAnchor.que}.{cnsAnchor.ref}")
        fm6_blastn_file = os.path.join(blastnDir, f"{cnsAnchor.que}.{cnsAnchor.ref}",
                                       f"{cnsAnchor.que}.{cnsAnchor.ref}.blastn.fmt6.txt")
        BasePan.pymkdir(outDir)
        track_progress = parallelScheduler.make_call_back()
        error_callback = parallelScheduler.make_error_callback()
        result = pool.apply_async(CNSsAnchor.multi_JCVI, args=(cnsAnchor, outDir, fm6_blastn_file,),
                                  callback=track_progress, error_callback=error_callback)
        results.append(result)
    pool.close()
    pool.join()
    del total_tasks, completed_tasks, progress_block, parallelScheduler, results
    logger.info("parallel jcvi align Finish")

    # 4. parallel jcvi halLiftover filter
    logger.info("parallel parallel jcvi halLiftover filter start ...")
    JCVI_filter_threads = args_dict["jcvi_parallel"]
    total_tasks = len(CNSsAnchorList)
    logger.info("total_tasks: {}".format(total_tasks))
    logger.info("threads: {}".format(JCVI_filter_threads))
    completed_tasks = multiprocessing.Manager().Value('i', 0)
    results = []
    pool = multiprocessing.Pool(processes=JCVI_filter_threads)
    progress_block = JCVI_filter_threads
    parallelScheduler = ParallelScheduler(logger, total_tasks, completed_tasks, progress_block)
    for cnsAnchor in CNSsAnchorList:
        outDir = os.path.join(JCVIDir, f"{cnsAnchor.que}.{cnsAnchor.ref}")
        track_progress = parallelScheduler.make_call_back()
        error_callback = parallelScheduler.make_error_callback()
        result = pool.apply_async(CNSsAnchor.multi_halLiftoverFilter_JCVI, args=(cnsAnchor, outDir, workDir,),
                                  callback=track_progress, error_callback=error_callback)
        results.append(result)
    pool.close()
    pool.join()
    del total_tasks, completed_tasks, progress_block, parallelScheduler, results
    logger.info("parallel parallel jcvi halLiftover filter Finish")

    # 5. parallel rescue cns anchor
    logger.info("parallel rescue cns anchor start ...")
    rescue_threads = args_dict["jcvi_parallel"]
    total_tasks = len(CNSsAnchorList)
    logger.info("total_tasks: {}".format(total_tasks))
    logger.info("threads: {}".format(rescue_threads))
    completed_tasks = multiprocessing.Manager().Value('i', 0)
    results = []
    pool = multiprocessing.Pool(processes=JCVI_filter_threads)
    progress_block = JCVI_filter_threads
    parallelScheduler = ParallelScheduler(logger, total_tasks, completed_tasks, progress_block)
    for cnsAnchor in CNSsAnchorList:
        track_progress = parallelScheduler.make_call_back()
        error_callback = parallelScheduler.make_error_callback()
        result = pool.apply_async(run_rescue_cns_anchor, args=(cnsAnchor, geneWorkDir, configData,),
                                  callback=track_progress, error_callback=error_callback)
        results.append(result)
    pool.close()
    pool.join()
    del total_tasks, completed_tasks, progress_block, parallelScheduler, results
    logger.info("parallel rescue cns anchor Finish")

    # 6. merge jcvi and rescue cns anchor
    logger.info("merge jcvi and rescue cns anchor start ...")
    for cnsAnchor in CNSsAnchorList:
        rescued_map_file = os.path.join(cnsAnchor.workDir, "JCVIDir", f"{cnsAnchor.que}.{cnsAnchor.ref}",
                                        f"{cnsAnchor.que}.{cnsAnchor.ref}.rescue_cns_map.csv")
        out_file_1 = os.path.join(cnsAnchor.workDir, "JCVIDir", f"{cnsAnchor.que}.{cnsAnchor.ref}",
                                  f"{cnsAnchor.que}.{cnsAnchor.ref}.halLiftoverFilter.rescue.anchors")
        out_file_2 = os.path.join(cnsAnchor.workDir, "JCVIDir", f"{cnsAnchor.que}.{cnsAnchor.ref}",
                                  f"{cnsAnchor.que}.{cnsAnchor.ref}.halLiftoverFilter.rescue.lifted.anchors")
        BasePan.merge_files(out_file_1, rescued_map_file,
                            os.path.join(cnsAnchor.workDir, "JCVIDir", f"{cnsAnchor.que}.{cnsAnchor.ref}",
                                         f"{cnsAnchor.que}.{cnsAnchor.ref}.halLiftoverFilter.anchors"))
        BasePan.merge_files(out_file_2, rescued_map_file,
                            os.path.join(cnsAnchor.workDir, "JCVIDir", f"{cnsAnchor.que}.{cnsAnchor.ref}",
                                         f"{cnsAnchor.que}.{cnsAnchor.ref}.halLiftoverFilter.lifted.anchors"))
    logger.info("merge jcvi and rescue cns anchor Finish !!!")

    # 7. parallel get Two-way anchor
    logger.info("parallel get Two-way anchor start ...")
    TwoWay_parallel_threads = args_dict["jcvi_parallel"]
    total_tasks = len(CNSsAnchorList)
    logger.info("total_tasks: {}".format(total_tasks))
    logger.info("threads: {}".format(TwoWay_parallel_threads))
    completed_tasks = multiprocessing.Manager().Value('i', 0)
    results = []
    pool = multiprocessing.Pool(processes=TwoWay_parallel_threads)
    progress_block = TwoWay_parallel_threads
    parallelScheduler = ParallelScheduler(logger, total_tasks, completed_tasks, progress_block)
    for cnsAnchor in CNSsAnchorList:
        track_progress = parallelScheduler.make_call_back()
        error_callback = parallelScheduler.make_error_callback()
        result = pool.apply_async(CNSsAnchor.multi_get_TowWey, args=(cnsAnchor,),
                                  callback=track_progress, error_callback=error_callback)
        results.append(result)
    pool.close()
    pool.join()
    del total_tasks, completed_tasks, progress_block, parallelScheduler, results
    logger.info("parallel get Two-way anchor Finish")

    logger.info("End mergeCNSsMap !!!")

