import os
import glob
import pandas as pd
from panCG.lib.base import BasePan
from panCG.lib.pancns import MergeSplit


def muli_MergeSplit(data_a, data_b, merge_col_on):
    MergeSpliter = MergeSplit(data_a, data_b, merge_col_on)
    result_series = MergeSpliter.merge_data()
    return result_series


def run_geneClustering(logger, config, workDir, Reference, args_dict):
    """
    Clustering using orthofinder
    Args:
        args_dict:
        logger:
        config:
        workDir:
        Reference:
    Returns:
        None
    """
    logger.info("--------------------------- step 2. Start orthofinder geneClustering ... -----------------------------")
    configData = BasePan.read_yaml(config)
    species_tree = configData["species_tree"]
    speciesList = BasePan.get_species_li(species_tree, Reference)
    OrthoFinderDir = os.path.join(workDir, "OrthoFinderDir")
    BasePan.pymkdir(OrthoFinderDir)

    orthofinder_ok_file = os.path.join(OrthoFinderDir, "orthofinder.ok")
    if not os.path.exists(orthofinder_ok_file):
        logger.info("OrthoFinderDir: {}".format(OrthoFinderDir))
        for i in speciesList:
            longest_pep_file = configData["species"][i]["longest_pep_fasta"]
            os.symlink(longest_pep_file, os.path.join(OrthoFinderDir, f"{i}.longest.pep.fasta"))

        logger.info("Run orthofinder ...")
        orthofinder_cmd = BasePan.check_cmd("orthofinder")
        OrthoFinderThreads = args_dict["orthofinder_threads"]
        cmd = f"{orthofinder_cmd} -f {OrthoFinderDir} -t {OrthoFinderThreads} > orthofinder.log 2>&1"
        BasePan.cmd_linux(cmd)
        logger.info("Run orthofinder Finish")
        BasePan.touch_file(orthofinder_ok_file)
    else:
        logger.info("jump orthofinder ...")
    del orthofinder_ok_file

    RefFreePanGene = glob.glob(os.path.join(OrthoFinderDir, "OrthoFinder", "*", "Orthogroups", "Orthogroups.tsv"),
                               recursive=True)[0]
    logger.info("groups data: {}".format(RefFreePanGene))

    RefFreePanGeneData = pd.read_csv(RefFreePanGene, sep="\t")
    RefFreePanGeneData = RefFreePanGeneData.fillna('.')
    RefFreePanGeneData.columns = RefFreePanGeneData.columns.map(lambda x: x.split('.')[0])
    newColumnOrder = ["Orthogroup"]
    newColumnOrder.extend(speciesList)
    RefFreePanGeneData = RefFreePanGeneData[newColumnOrder]
    # The results of OrthoFinder use ", " to separate genes.
    # Now replace ", " with ", "
    RefFreePanGeneData = RefFreePanGeneData.replace({', ': ','}, regex=True)
    logger.info("speciesList: {}".format(RefFreePanGeneData.columns.tolist()))

    # Add the genes that are not clustered by OrthoFinder
    logger.info("Add the genes that are not clustered by OrthoFinder")
    unMapGeneData_li = []
    for specie in speciesList:
        logger.info(f"add {specie} ...")
        allGeneFile = configData["species"][specie]["longest_pep_bed"]
        bedData = pd.read_table(allGeneFile, sep='\t', names=['chrID', 'start', 'end', 'geneID', 'score', 'chain'])
        allGeneList = bedData["geneID"].to_list()
        mapGeneList = [i.strip() for item in RefFreePanGeneData[specie].to_list() if item != "." for i in item.split(",")]
        unMapGeneList = list(set(allGeneList) - set(mapGeneList))
        unMapGeneData_li.append(pd.DataFrame({specie: unMapGeneList}))
    unMapGeneData = pd.concat(unMapGeneData_li, axis=0).fillna(".")
    new_column = ["UnMapOG{:0{}d}".format(i, len(str(len(unMapGeneData))) + 1) for i in
                  range(1, len(unMapGeneData) + 1)]
    unMapGeneData.insert(0, "Orthogroup", new_column)
    RefFreePanGeneData = pd.concat([RefFreePanGeneData, unMapGeneData], axis=0, ignore_index=True)
    RefFreePanGeneData = RefFreePanGeneData.rename(columns={'Orthogroup': 'Group'})
    ClusterDir = os.path.join(workDir, "Cluster")
    BasePan.pymkdir(ClusterDir)
    indexFile = os.path.join(ClusterDir, "All.Cluster.csv")
    RefFreePanGeneData.to_csv(indexFile, sep='\t', index=False)
    logger.info("geneClustering End")


