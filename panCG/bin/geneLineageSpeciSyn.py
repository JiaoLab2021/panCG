import pandas as pd
import os
import glob
import networkx as nx
import multiprocessing

from panCG.lib.base import ParallelScheduler
from panCG.lib.base import BasePan
from panCG.lib.PhySynNet import BulitPhySynNet
from panCG.lib.LineageSpecicSynteny import GeneLineageSpecicSynteny
from panCG.lib.LineageSpecicSynteny import GeneLineageSpeSynNoOrth


def glss(logger, config, workDir, Reference, lineage_species_file, output_file, args_dict, threads):
    logger.info(f"Parameter: {args_dict}")

    no_orthofinder = args_dict["no_orthofinder"]  # default value of True
    configData = BasePan.read_yaml(config)
    species_tree = configData["species_tree"]
    speciesList = BasePan.get_species_li(species_tree, Reference)

    df = pd.read_csv(lineage_species_file, header=None, sep="\t")
    Lineage_species = list(set(df[0].tolist()))
    del df
    logger.info("Lineage species: {}".format(",".join(Lineage_species)))

    logger.info("Reading Synteny gene graph network ...")
    GeneLineageSpecicSynteny.set_k(args_dict["k_clique"])
    GeneLineageSpecicSynteny.set_miss_spe_num(args_dict["miss_spe_num"])
    graph_file = os.path.join(workDir, "JCVIDir", "Synteny.gene.overlap.graph.pkl")
    BulitPhySynNeter = BulitPhySynNet(workDir, graph_file, speciesList)
    G = BulitPhySynNeter.read_Synteny_Network()

    if no_orthofinder:
        logger.info("Reading OrthoFinder Orthogroups ...")
        OrthoFinderDir = os.path.join(workDir, "OrthoFinderDir")
        Orthogroups_file = glob.glob(os.path.join(OrthoFinderDir, "OrthoFinder", "*", "Orthogroups", "Orthogroups.tsv"), recursive=True)[0]
        Orthogroups_data = pd.read_csv(Orthogroups_file, sep="\t")
        Orthogroups_data = Orthogroups_data.fillna('.')
        Orthogroups_data.columns = Orthogroups_data.columns.map(lambda x: x.split('.')[0])
        newColumnOrder = ["Orthogroup"]
        newColumnOrder.extend(speciesList)
        Orthogroups_data = Orthogroups_data[newColumnOrder]
        Orthogroups_data = Orthogroups_data.replace({', ': ','}, regex=True)
        Orthogroups_data.set_index("Orthogroup", inplace=True)

        total_tasks = threads * 5
        chunk_num = len(Orthogroups_data) // total_tasks
        logger.info("total_tasks: {}".format(total_tasks))
        logger.info("threads: {}".format(threads))
        completed_tasks = multiprocessing.Manager().Value('i', 0)
        results = []
        pool = multiprocessing.Pool(processes=threads)
        progress_block = threads
        parallelScheduler = ParallelScheduler(logger, total_tasks, completed_tasks, progress_block)

        for i in range(0, len(Orthogroups_data), chunk_num):
            track_progress = parallelScheduler.make_call_back()
            error_callback = parallelScheduler.make_error_callback()

            chunk_data = Orthogroups_data.iloc[i:i + chunk_num]
            LineageSpecicSyntenyer = GeneLineageSpecicSynteny(G, Lineage_species, chunk_data)
            result = pool.apply_async(LineageSpecicSyntenyer.run,
                                      callback=track_progress, error_callback=error_callback)
            results.append(result)
        pool.close()
        pool.join()

        data = pd.concat([i.get()[0] for i in results], axis=0, ignore_index=True)
        all_data = pd.concat([i.get()[1] for i in results], axis=0, ignore_index=True)
        all_data.to_csv("orthofinder.cluster.cpm.tsv", sep="\t", index=False)
    else:
        # Get the species corresponding to each gene
        gene2dict = {}
        for spe in speciesList:
            gene_bed_file = configData["species"][spe]["longest_pep_bed"]
            gene_bed_df = pd.read_csv(gene_bed_file, sep="\t", header=None)
            dict_tmp = {i: spe for i in gene_bed_df[3].tolist()}
            gene2dict.update(dict_tmp)
        connected_components = list(nx.connected_components(G))

        total_tasks = threads * 5
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

            set_li = connected_components[i:i + chunk_num]
            merged_set = set().union(*set_li)
            sub_G_ = G.subgraph(merged_set).copy()

            GeneLineageSpeSynNoOrther = GeneLineageSpeSynNoOrth(sub_G_, Lineage_species)
            result = pool.apply_async(GeneLineageSpeSynNoOrther.run, args=(gene2dict,),
                                      callback=track_progress, error_callback=error_callback)
            results.append(result)
        pool.close()
        pool.join()
        all_cluster_data = pd.concat([i.get()[0] for i in results], axis=0, ignore_index=True)
        all_cluster_data.to_csv("no_orthofinder.cluster.cpm.tsv", sep="\t", index=False)
        data = pd.concat([i.get()[1] for i in results], axis=0, ignore_index=True)

    data.to_csv(output_file, sep="\t", index=False)


