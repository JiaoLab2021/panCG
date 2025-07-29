import os
import pandas as pd
import multiprocessing
from panCG.lib.base import BasePan
from panCG.lib.base import ParallelScheduler
from panCG.lib.sqliteDatabase import SQLiteDatabase
from panCG.lib.indexAssign import GeneIndex
from panCG.lib.CPM import CPM
from panCG.lib.PhySynNet import BulitPhySynNet


def parallel_GeneIndex(diamond_blastp_db, chunk, G, k):
    GeneIndexer = GeneIndex(diamond_blastp_db, chunk, G, k)
    return GeneIndexer.run()


def run_pangeneAssign(logger, config, workDir, Reference, args_dict):
    logger.info("--------------------------------- step 3. Start geneIndexAssign ... --------------------------------")
    configData = BasePan.read_yaml(config)
    species_tree = configData["species_tree"]
    speciesList = BasePan.get_species_li(species_tree, Reference)
    diamondDir = os.path.join(workDir, "diamondDir")

    max_cliques_num = args_dict["max_cliques_num"]
    CPM.set_max_cliques_num(max_cliques_num)

    # diamond_blastp data import database
    logger.info("diamond_blastp data import database ...")
    diamond_blastp_db_name = os.path.join(workDir, 'diamondDir', 'diamond_blastp.db')
    index_column = 'que'  # index col
    data_columns = ['ref', 'score']  # data col
    diamond_blastp_db = SQLiteDatabase(diamond_blastp_db_name, index_column, data_columns)
    if not os.path.exists(diamond_blastp_db_name):
        diamond_blastp_df_li = []
        for i in range(len(speciesList)):
            for j in range(len(speciesList)):
                if i == j:
                    continue
                que, ref = speciesList[i], speciesList[j]
                diamond_blastp_fmt6_file = os.path.join(diamondDir, f"{que}.{ref}",
                                                        f"{que}.{ref}.diamond_blastp.filter.fmt6.txt")
                df = pd.read_csv(diamond_blastp_fmt6_file, sep="\t",
                                 header=None, usecols=[0, 1, 11], names=["que", "ref", "score"])
                diamond_blastp_df_li.append(df)
        diamond_blastp_db.insert_data(diamond_blastp_df_li)
        diamond_blastp_db.create_index()
        del diamond_blastp_df_li
    else:
        pass
    del index_column, data_columns

    cluster_file = os.path.join(workDir, "Cluster", "All.Cluster.csv")
    Orthogroups_data = pd.read_csv(cluster_file, sep="\t")
    Orthogroups_data = Orthogroups_data.set_index("Group")
    Orthogroups_data = Orthogroups_data[speciesList]

    # reading Phylogenomic Synteny Network
    logger.info("reading Phylogenomic Synteny Network ...")
    graph_file = os.path.join(workDir, "JCVIDir", "Synteny.gene.overlap.graph.pkl")
    BulitPhySynNeter = BulitPhySynNet(workDir, graph_file, speciesList)
    G = BulitPhySynNeter.read_Synteny_Network()

    threads = args_dict["assign_threads"]
    chunk_size = args_dict["assign_chunk_size"]
    k_clique = args_dict["k_clique"]
    GeneIndex.set_single_abstract_cutoff(args_dict["single_abstract_cutoff"])

    chunks = [Orthogroups_data.iloc[i:i + chunk_size] for i in range(0, len(Orthogroups_data), chunk_size)]
    total_tasks = len(chunks)
    logger.info("total_tasks: {}".format(total_tasks))
    logger.info("threads: {}".format(threads))
    completed_tasks = multiprocessing.Manager().Value('i', 0)
    results = []
    pool = multiprocessing.Pool(processes=threads)
    progress_block = threads
    parallelScheduler = ParallelScheduler(logger, total_tasks, completed_tasks, progress_block)
    for chunk in chunks:
        track_progress = parallelScheduler.make_call_back()
        error_callback = parallelScheduler.make_error_callback()
        result = pool.apply_async(parallel_GeneIndex, args=(diamond_blastp_db, chunk, G, k_clique,),
                                  callback=track_progress, error_callback=error_callback)
        results.append(result)
    pool.close()
    pool.join()

    df_li = []
    for i in results:
        df_li.append(i.get())
    result = pd.concat(df_li, axis=0, ignore_index=True)
    result = result.fillna('.')

    out_col_li = ["Group", "Index"]
    out_col_li.extend(speciesList)
    result = result[out_col_li]

    RefIndexDir = os.path.join(workDir, "Ref_{}_IndexDir".format(Reference))
    BasePan.pymkdir(RefIndexDir)
    gene_index_assign_file = os.path.join(RefIndexDir, "Ref.{}.cpm.cluster.csv".format(Reference))
    result.to_csv(gene_index_assign_file, sep="\t", index=False)
    logger.info("geneIndexAssign End")

