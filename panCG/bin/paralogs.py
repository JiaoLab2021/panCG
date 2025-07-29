import os
import glob
import pandas as pd
from panCG.lib.base import BasePan
from panCG.lib.pangeneIndex import ParseTree


def parse_paralogs(logger, config, workDir, Reference, args_dict):
    logger.info("--------------------------------- step 4. Start parse paralogs ... --------------------------------")

    RefIndexDir = os.path.join(workDir, "Ref_{}_IndexDir".format(Reference))
    gene_index_assign_file = os.path.join(RefIndexDir, "Ref.{}.cpm.cluster.csv".format(Reference))
    geneIndexMergeData = pd.read_csv(gene_index_assign_file, sep="\t")

    logger.info("parse paralogs ...")
    min_paralogs_spe_num = args_dict["min_paralogs_spe_num"]
    ParseTree.set_cutoff(min_paralogs_spe_num)
    configData = BasePan.read_yaml(config)
    species_tree = configData["species_tree"]
    speciesList = BasePan.get_species_li(species_tree, Reference)

    OrthoFinderDir = os.path.join(workDir, "OrthoFinderDir")
    tree_dir = glob.glob(os.path.join(OrthoFinderDir, "OrthoFinder", "*", "Gene_Trees"), recursive=True)[0]
    data_li = []
    all_len = len(geneIndexMergeData)
    for index, row in geneIndexMergeData.iterrows():
        if index % 1000 == 0:
            schedule_ = (index / all_len) * 100
            logger.info("{:.2f} % completed".format(schedule_))
        name = row["Group"]
        tree_file = os.path.join(tree_dir, f"{name}_tree.txt")
        if os.path.exists(tree_file):
            treer = ParseTree(row, speciesList, tree_file)
            df = treer.run()
            data_li.append(df)
        else:
            data_li.append(pd.DataFrame(row).transpose())
    result = pd.concat(data_li, axis=0, ignore_index=True)
    out_col_li = ["Group", "Index"]
    out_col_li.extend(speciesList)
    result = result[out_col_li]
    result.to_csv(os.path.join(RefIndexDir, "Ref.{}.panGene.final.csv".format(Reference)), sep='\t', header=True, index=False)
    logger.info("parse paralogs End !!!")


