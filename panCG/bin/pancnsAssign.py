import pandas as pd
import os
import multiprocessing

from panCG.lib.base import BasePan
from panCG.lib.base import ParallelScheduler
from panCG.lib.pancns import read_map
from panCG.lib.indexAssign import multi_IndexAssign


def run_pancnsAssign(logger, config, workDir, Reference, args_dict):
    logger.info("--------------------------------- step 3. Start cnsIndexAssign ... ----------------------------------")
    configData = BasePan.read_yaml(config)
    species_tree = configData["species_tree"]
    speciesList = BasePan.get_species_li(species_tree, Reference)

    JCVIDir = os.path.join(workDir, "JCVIDir")

    # 2. Store all synteny CNS in the dictionary
    logger.info('Store all synteny CNS in the dictionary ...')
    JCVI_Dict = {}
    for que in speciesList:
        for ref in speciesList:
            if que == ref:
                continue
            jcvi_filter_File = os.path.join(JCVIDir, f"{que}.{ref}",
                                            f"{que}.{ref}.halLiftoverFilter.rescue.lifted.anchors")
            JCVI_Dict[f"{que}.{ref}"] = read_map(jcvi_filter_File)

    logger.info("Finish data loading")

    cnsCluster_file = os.path.join(workDir, "Index", "cnsCluster.csv")
    df = pd.read_csv(cnsCluster_file, sep="\t")
    random_seed = args_dict["assign_random_seed"]
    # Disorganize the data frame and add random seeds for repeatability
    df = df.sample(frac=1, random_state=random_seed).reset_index(drop=True)
    speciesList = df.columns.tolist()[1:]
    # new_column = ["Group{:0{}d}".format(i, len(str(len(df))) + 1) for i in range(1, len(df) + 1)]
    # df.insert(0, "Group", new_column)

    chunk_size = args_dict["assign_chunk_size"]
    threads = args_dict["assign_threads"]

    total_tasks = len(range(0, len(df), chunk_size))
    logger.info("total_tasks: {}".format(total_tasks))  
    logger.info("threads: {}".format(threads))  
    completed_tasks = multiprocessing.Manager().Value('i', 0)
    results = []
    pool = multiprocessing.Pool(processes=threads)
    progress_block = threads
    parallelScheduler = ParallelScheduler(logger, total_tasks, completed_tasks, progress_block)
    for i in range(0, len(df), chunk_size):
        track_progress = parallelScheduler.make_call_back()
        error_callback = parallelScheduler.make_error_callback()
        chunk_data = df.iloc[i:i + chunk_size]
        result = pool.apply_async(multi_IndexAssign, args=(chunk_data, speciesList, JCVI_Dict,), callback=track_progress, error_callback=error_callback)
        results.append(result)
    pool.close()
    pool.join()

    Ref_Speci_Data = pd.concat([i.get() for i in results], axis=0)

    RefIndexDir = os.path.join(workDir, "Ref_{}_IndexDir".format(Reference))
    BasePan.pymkdir(RefIndexDir)

    cnsIndexAssignFile = os.path.join(RefIndexDir, "Ref.{}.cnsIndexAssign.csv".format(Reference))
    Ref_Speci_Data.to_csv(cnsIndexAssignFile, sep='\t', index=False)

    logger.info("End cnsIndexAssign !!!")


