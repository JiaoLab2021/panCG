import pandas as pd
import os
import copy
import networkx as nx
import pickle

from panCG.lib.base import BasePan
from panCG.bin.cnsMapMerge import run_cnsMapMerge
from panCG.lib.pancnsFilter import PanCnsFilter
from panCG.lib.LineageSpecicSynteny import CnsPhySynNet


def run_CNS_synteny_net(logger, panCns, cnsConfig, workDir, reference, geneWorkDir, args_dict):
    logger.info("----------------------------------- CNS synteny network ... ------------------------------------")
    logger.info("args: {}".format(args_dict))
    cnsConfigData = BasePan.read_yaml(cnsConfig)
    species_tree = cnsConfigData["species_tree"]
    speciesList = BasePan.get_species_li(species_tree, reference)
    logger.info("species list: {}".format(", ".join(speciesList)))

    BasePan.pymkdir(workDir)
    os.chdir(workDir)

    new_cnsConfigData = copy.deepcopy(cnsConfigData)
    new_cnsConfig_file = os.path.join(workDir, "new_cns.yaml")
    cns_dir = os.path.join(workDir, "cns_dir")
    BasePan.pymkdir(cns_dir)
    panCns_data = pd.read_csv(panCns, sep="\t")
    panCnsFilter = PanCnsFilter()
    for spe in speciesList:
        logger.info(f"\tparse {spe} ...")
        all_cns = panCnsFilter.ToBed(panCns_data, "Index", spe)
        out_cns_file = os.path.join(cns_dir, f"{spe}.cns_index.bed")
        all_cns.to_csv(out_cns_file, sep="\t", index=False, header=False)
        new_cnsConfigData["species"][spe]["CNS_bed"] = out_cns_file
    BasePan.write_yaml(new_cnsConfigData, new_cnsConfig_file)

    logger.info(f"Run CNS Map ...")
    run_cnsMapMerge(logger, new_cnsConfig_file, workDir, geneWorkDir, reference, args_dict)
    logger.info(f"Finish CNS Map ")

    logger.info(f"Run build CNS PhySynNet ...")
    graph_file = os.path.join(workDir, "JCVIDir", "Synteny.CNS.graph.pkl")
    CnsPhySynNeter = CnsPhySynNet(workDir, graph_file, speciesList)
    CnsPhySynNeter.build_net()

    # logger.info(f"loading CNS Map result ...")
    # df_li = []
    # for i in range(len(speciesList)):
    #     logger.info(f"\tparse {i} ...")
    #     for j in range(i + 1, len(speciesList)):
    #         que, ref = speciesList[i], speciesList[j]
    #         if que == ref:
    #             continue
    #         net_file = os.path.join(workDir, "JCVIDir", f"{que}.{ref}", f"{que}.{ref}.halLiftoverFilter.rescue.lifted.anchors")
    #         df = pd.read_csv(net_file, sep="\t", header=None, names=["source", "target"])
    #         df["source"] = df["source"].apply(lambda x: f"{que}:{x}")
    #         df["target"] = df["target"].apply(lambda x: f"{ref}:{x}")
    #         df_li.append(df)
    # data = pd.concat(df_li, axis=0, ignore_index=True)
    # del df_li
    #
    # G = nx.from_pandas_edgelist(data, source='source', target='target')
    # graph_file = os.path.join(workDir, "JCVIDir", "Synteny.CNS.graph.pkl")
    # with open(graph_file, "wb") as fo:
    #     pickle.dump(G, fo)

    logger.info(f"The Synteny CNS network is located at: {graph_file}")

