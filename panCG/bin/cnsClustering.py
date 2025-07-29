import pandas as pd
import os
import pickle
import networkx as nx
from panCG.lib.base import BasePan


def run_cnsClustering(logger, config, workDir, reference):
    """
    This function uses the results of halLiftover filtering to cluster cns
    Args:
        reference:
        logger: logger
        config: config
        workDir: workDir
    Returns:
        None
    """
    logger.info("--------------------------------- step 2. Start cnsClustering ... -----------------------------------")
    configData = BasePan.read_yaml(config)
    IndexDir = os.path.join(workDir, "Index")
    BasePan.pymkdir(IndexDir)
    halLiftoverDir = os.path.join(workDir, "halLiftoverDir")
    species_tree = configData["species_tree"]
    speciesList = BasePan.get_species_li(species_tree, reference)
    graph_file = os.path.join(IndexDir, "halLiftover_anchors.net.pickle")
    logger.info("Generate clustering network ...")
    G = nx.Graph()
    for que in speciesList:
        logger.info(f"\tProcessing {que} ...")
        for ref in speciesList:
            if que == ref:
                continue
            anchorFile = os.path.join(halLiftoverDir, f"{que}.{ref}", f"{que}.{ref}.halLiftover.anchors")
            anchorData = pd.read_csv(anchorFile, sep='\t', comment='#', usecols=[0, 1], names=[que, ref])
            anchorData[que] = anchorData[que].apply(lambda x: f"{que}::{x}")
            anchorData[ref] = anchorData[ref].apply(lambda x: f"{ref}::{x}")
            graph = nx.from_pandas_edgelist(anchorData, source=que, target=ref)
            G.update(graph)
    with open(graph_file, "wb") as fo:
        pickle.dump(G, fo)
    logger.info(f"The clustering network is saved in {graph_file}")

    logger.info("Clustering based on connectivity graph ...")
    connected_components = list(nx.connected_components(G))
    logger.info("Found {} clusters".format(len(connected_components)))
    data_li = []
    for components in connected_components:
        dict_ = {}
        for i in components:
            spe, cns = i.split("::")[0], i.split("::")[1]
            dict_.setdefault(spe, set()).add(cns)
        dict_ = {k: ",".join(v) for k, v in dict_.items()}
        data_li.append(dict_)
    data = pd.DataFrame(data_li)
    data = data.fillna(".")

    logger.info("Add unclustered cns ...")
    df_li = [data]
    for specie in speciesList:
        cns_file = configData["species"][specie]["CNS_bed"]
        cns_data = pd.read_csv(cns_file, sep="\t", header=None, names=["chrID", "start", "end", "cnsID"])
        all_cns_set = set(cns_data["cnsID"].tolist())
        net_cns_set = {
            item
            for element in data[specie].tolist() if element != "."
            for item in element.split(",")
        }
        other_cns_set = all_cns_set - net_cns_set
        logger.info(f"\t{specie}\t{len(other_cns_set)}")
        df = pd.DataFrame(other_cns_set, columns=[specie])
        df_li.append(df)
    result = pd.concat(df_li, axis=0, ignore_index=True)
    result = result.fillna(".")
    result = result[speciesList]

    new_column = ["Group{:0{}d}".format(i, len(str(len(result))) + 1) for i in range(1, len(result) + 1)]
    result.insert(0, "Group", new_column)

    indexFile = os.path.join(IndexDir, "cnsCluster.csv")
    result.to_csv(indexFile, sep="\t", index=False)

    logger.info("End cnsClustering !!!")

