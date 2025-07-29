import pandas as pd
import pickle
import networkx as nx
import multiprocessing

from panCG.lib.base import ParallelScheduler
from panCG.lib.LineageSpecicSynteny import CnsLineageSpecicSynteny


def clss(logger, net_file, lineage_species_file, output_file, args_dict, threads):
    logger.info("args: {}".format(args_dict))

    df = pd.read_csv(lineage_species_file, header=None, sep="\t")
    Lineage_species = list(set(df[0].tolist()))
    del df
    logger.info("Lineage species: {}".format(",".join(Lineage_species)))

    CnsLineageSpecicSynteny.set_k(args_dict["k_clique"])
    CnsLineageSpecicSynteny.set_miss_spe_num(args_dict["miss_spe_num"])

    logger.info("loading CNS PhySynNet ...")
    with open(net_file, "rb") as f:
        G = pickle.load(f)

    logger.info("Clustering based on connectivity graph ...")
    connected_components = list(nx.connected_components(G))
    logger.info("Found {} clusters".format(len(connected_components)))
    components_index_li = ["SynGroup{:0{}d}".format(i, len(str(len(connected_components))) + 1) for i in range(1, len(connected_components) + 1)]

    data_li = []
    for components_index, components in zip(components_index_li, connected_components):
        dict_ = {}
        for i in components:
            spe, cns = i.split("::")[0], i.split("::")[1]
            dict_.setdefault(spe, set()).add(cns)
        dict1_ = {k: ",".join(v) for k, v in dict_.items()}
        dict1_["SynGroup"] = components_index
        data_li.append(dict1_)
    cluster_data = pd.DataFrame(data_li)
    cluster_data = cluster_data.fillna(".")
    cluster_data.set_index('SynGroup', inplace=True)
    cluster_data.to_csv("connected_components.data.tsv", sep="\t", index=False)
    print(cluster_data)
    logger.info(f"Generate cluster data， data.shape: {cluster_data.shape}")

    total_tasks = threads * 30
    chunk_num = len(connected_components) // total_tasks
    logger.info("total_tasks: {}".format(total_tasks))
    logger.info("threads: {}".format(threads))
    completed_tasks = multiprocessing.Manager().Value('i', 0)
    results = []
    pool = multiprocessing.Pool(processes=threads)
    progress_block = threads
    parallelScheduler = ParallelScheduler(logger, total_tasks, completed_tasks, progress_block)

    for i in range(0, len(connected_components), chunk_num):
        track_progress = parallelScheduler.make_call_back()
        error_callback = parallelScheduler.make_error_callback()

        chunk_index_li = components_index_li[i: i + chunk_num]
        chunk_cluster_li = connected_components[i: i + chunk_num]
        G_sub = G.subgraph(set.union(*chunk_cluster_li)).copy()

        chunk_cluster_dict = {}
        for components_index, components in zip(chunk_index_li, chunk_cluster_li):
            chunk_cluster_dict[components_index] = components
        CnsLineageSpecicSyntenyer = CnsLineageSpecicSynteny(G_sub, Lineage_species, chunk_cluster_dict)
        result = pool.apply_async(CnsLineageSpecicSyntenyer.run,
                                  callback=track_progress, error_callback=error_callback)
        results.append(result)
    pool.close()
    pool.join()

    data = pd.concat([i.get()[0] for i in results], axis=0, ignore_index=True)
    all_data = pd.concat([i.get()[1] for i in results], axis=0, ignore_index=True)

    logger.info("The size of Lineage Specic data is {}".format(data.shape))
    logger.info("The size of all cluster data is {}".format(all_data.shape))

    data.to_csv(output_file, sep="\t", index=False)
    all_data.to_csv("all_cns.cluster.tsv", sep="\t", index=False)





