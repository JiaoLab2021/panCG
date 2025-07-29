import os
import io
import pandas as pd
from panCG.lib.base import BasePan
import copy
from ete3 import Tree


class NoDirException(Exception):
    """ Custom exception types """
    def __init__(self, message):
        super().__init__(message)


class GeneAnchor(BasePan):
    diamond_run_threads = 6
    jcvi_run_threads = 6
    RBH_run_threads = 6
    JCVI_ENV = {}
    identity_threshold = 50  # %
    E_value = 1e-10

    diamond_cmd = BasePan.check_cmd("diamond")

    def __init__(self, que, ref, queProt, refProt, queLongestProtBed, refLongestProtBed, workDir):
        self.que, self.ref = que, ref
        self.queProt, self.refProt = queProt, refProt
        self.queLongestProtBed, self.refLongestProtBed = queLongestProtBed, refLongestProtBed
        self.workDir = workDir
        self.checking()  # check LongestProtDir and JCVIDir directory

    @classmethod
    def set_diamond_run_threads(cls, diamond_run_threads_num):
        cls.diamond_run_threads = diamond_run_threads_num

    @classmethod
    def set_jcvi_run_threads(cls, jcvi_run_threads_num):
        cls.jcvi_run_threads = jcvi_run_threads_num

    @classmethod
    def set_E_value(cls, Evalue):
        cls.E_value = Evalue

    @classmethod
    def set_identity_threshold(cls, identity_threshold_value):
        cls.identity_threshold = identity_threshold_value

    def checking(self):
        Dir = os.path.join(self.workDir, "JCVIDir")
        if not os.path.exists(Dir) or not os.path.isdir(Dir):
            raise NoDirException(f"The directory '{Dir}' does not exist.")
        Dir = os.path.join(self.workDir, "diamondDir")
        if not os.path.exists(Dir) or not os.path.isdir(Dir):
            raise NoDirException(f"The directory '{Dir}' does not exist.")
        Dir = os.path.join(self.workDir, "Merge_JCVI_diamond")
        if not os.path.exists(Dir) or not os.path.isdir(Dir):
            raise NoDirException(f"The directory '{Dir}' does not exist.")

    def Run_diamond(self):
        cmd = (f"{self.diamond_cmd} makedb --in {self.refProt} -d {self.ref}.diamond_db "
               f"> {self.que}.{self.ref}.diamond_makedb.log 2>&1")
        BasePan.cmd_linux(cmd)
        cmd = (f"{self.diamond_cmd} blastp -q {self.queProt} -d {self.ref}.diamond_db "
               f"-o {self.que}.{self.ref}.diamond_blastp.fmt6.txt "
               f"--very-sensitive --threads {self.diamond_run_threads} "
               f"> {self.que}.{self.ref}.diamond_blastp.log 2>&1")
        BasePan.cmd_linux(cmd)

    def Filter_diamond(self):
        """ Filter by identity """
        rawFile = f"{self.que}.{self.ref}.diamond_blastp.fmt6.txt"
        outFile = f"{self.que}.{self.ref}.diamond_blastp.filter.fmt6.txt"
        data = pd.read_csv(rawFile, sep="\t", header=None, comment='#', skip_blank_lines=True)
        data.columns = ["Query_id", "Target_id", "identity", "Length", "Mismatches", "Gap_openings", "Query_start",
                        "Query_end", "Target_start", "Target_end", "E_value", "Bit_score"]
        result = data[(data['identity'] >= self.identity_threshold) & (data['E_value'] < self.E_value)]
        result.to_csv(outFile, sep="\t", header=False, index=False)

    def Run_JCVI(self):
        if not os.path.exists(os.path.join(self.workDir, "JCVIDir", f"{self.que}.{self.ref}", f"{self.que}.bed")):
            os.symlink(
                self.queLongestProtBed,
                os.path.join(self.workDir, "JCVIDir", f"{self.que}.{self.ref}", f"{self.que}.bed")
            )
        if not os.path.exists(os.path.join(self.workDir, "JCVIDir", f"{self.que}.{self.ref}", f"{self.ref}.bed")):
            os.symlink(
                self.refLongestProtBed,
                os.path.join(self.workDir, "JCVIDir", f"{self.que}.{self.ref}", f"{self.ref}.bed")
            )
        if not os.path.exists(os.path.join(self.workDir, "JCVIDir", f"{self.que}.{self.ref}", f"{self.que}.pep")):
            os.symlink(
                self.queProt,
                os.path.join(self.workDir, "JCVIDir", f"{self.que}.{self.ref}", f"{self.que}.pep")
            )
        if not os.path.exists(os.path.join(self.workDir, "JCVIDir", f"{self.que}.{self.ref}", f"{self.ref}.pep")):
            os.symlink(
                self.refProt,
                os.path.join(self.workDir, "JCVIDir", f"{self.que}.{self.ref}", f"{self.ref}.pep")
            )
        if not os.path.exists(os.path.join(self.workDir, "JCVIDir", f"{self.que}.{self.ref}", f"{self.que}.{self.ref}.last")):
            os.symlink(
                os.path.join(self.workDir, "diamondDir", f"{self.que}.{self.ref}",
                            f"{self.que}.{self.ref}.diamond_blastp.filter.fmt6.txt"),
                os.path.join(self.workDir, "JCVIDir", f"{self.que}.{self.ref}", f"{self.que}.{self.ref}.last")
            )

        cmd = (f"python -m jcvi.compara.catalog ortholog --dbtype prot --min_size 5 --no_strip_names {self.que} {self.ref} "
               f"--cpus={self.jcvi_run_threads} > jcvi.compara.catalog.ortholog.log 2>&1")
        BasePan.cmd_linux(cmd)
        cmd = (f"python -m jcvi.compara.synteny screen --simple {self.que}.{self.ref}.anchors "
               f"{self.que}.{self.ref}.new > jcvi.compara.synteny.screen.log 2>&1")
        BasePan.cmd_linux(cmd)

    @staticmethod
    def parse_TowWey_fm6(in_file, out_file):
        data = pd.read_csv(in_file, sep="\t", header=None, comment='#', skip_blank_lines=True)
        columns = list(data.columns)
        columns[0], columns[1] = columns[1], columns[0]
        data = data[columns]
        data.to_csv(out_file, sep="\t", header=False, index=False)

    def get_TowWey_diamond(self):
        filter_raw_file = os.path.join(self.workDir, "diamondDir", f"{self.que}.{self.ref}",
                                       f"{self.que}.{self.ref}.diamond_blastp.filter.fmt6.txt")
        outDir = os.path.join(self.workDir, "diamondDir", f"{self.ref}.{self.que}")
        self.pymkdir(outDir)
        out_raw_file = os.path.join(outDir, f"{self.ref}.{self.que}.diamond_blastp.filter.fmt6.txt")
        self.parse_TowWey_fm6(filter_raw_file, out_raw_file)

    def parse_TowWey_jcvi(self, in_file, out_file):
        sequenceIO = io.StringIO()
        # txt = ""
        with open(in_file, "r") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    sequenceIO.write(line)
                else:
                    line_li = line.strip().split("\t")
                    txt = "{}\t{}\t{}\n".format(line_li[1], line_li[0], line_li[2])
                    sequenceIO.write(txt)
        result = sequenceIO.getvalue()
        sequenceIO.close()
        self.write_to(out_file, result)

    def get_TowWey_jcvi(self):
        anchors_file = os.path.join(self.workDir, "JCVIDir", f"{self.que}.{self.ref}", f"{self.que}.{self.ref}.anchors")
        lifted_anchors_file = os.path.join(self.workDir, "JCVIDir", f"{self.que}.{self.ref}", f"{self.que}.{self.ref}.lifted.anchors")
        outDir = os.path.join(self.workDir, "JCVIDir", f"{self.ref}.{self.que}")
        self.pymkdir(outDir)
        # anchors_file
        out_anchors_file = os.path.join(outDir, f"{self.ref}.{self.que}.anchors")
        self.parse_TowWey_jcvi(anchors_file, out_anchors_file)
        # lifted_anchors_file
        out_lifted_anchors_file = os.path.join(outDir, f"{self.ref}.{self.que}.lifted.anchors")
        self.parse_TowWey_jcvi(lifted_anchors_file, out_lifted_anchors_file)

    @staticmethod
    def get_anchor(in_file, out_file):
        try:
            data = pd.read_csv(in_file, sep="\t", header=None, comment='#', skip_blank_lines=True)
        except (pd.errors.EmptyDataError, FileNotFoundError):
            data = pd.DataFrame(columns=[0, 1])
        # data = pd.read_csv(in_file, sep="\t", header=None, comment='#', skip_blank_lines=True)
        data = data[[0, 1]]
        data.to_csv(out_file, sep="\t", header=False, index=False)

    @staticmethod
    def multi_RunDiamond(GeneAnchor_Instances, workDir):
        # Parallel interface
        os.chdir(workDir)
        GeneAnchor_Instances.Run_diamond()

    @staticmethod
    def multi_FilterDiamond(GeneAnchor_Instances, workDir):
        os.chdir(workDir)
        GeneAnchor_Instances.Filter_diamond()

    @staticmethod
    def multi_RunJCVI(GeneAnchor_Instances, workDir):
        # Parallel interface
        os.chdir(workDir)
        GeneAnchor_Instances.Run_JCVI()

    @staticmethod
    def multi_get_TowWey(GeneAnchor_Instances):
        GeneAnchor_Instances.get_TowWey_diamond()
        GeneAnchor_Instances.get_TowWey_jcvi()


class MulGeneMap:
    def __init__(self, anchorFile, que, ref, queGeneBed, refGeneBed, IndexDir, Data=pd.DataFrame(), CNSIndexData=None):
        self.Data = Data
        self.anchorFile = anchorFile
        self.que, self.ref = que, ref
        self.queGeneBed, self.refGeneBed = queGeneBed, refGeneBed
        self.IndexDir = IndexDir
        self.CNSIndexData = CNSIndexData

    def cluster_li(self, col_name):
        """
        This function aggregates the row indices of cells in a column of cells in a data frame that have repeated
        elements (separated by ,).
        :param col_name: Column Name
        :return: dict: {'1': ['Charlie'], '0-2': ['Alice'], '3': ['Charlie1']} and a list of row indices with null values
        """
        clusterIndexDir = {}
        emptyIndexLi = []
        cellLi = [i.split(",") for i in self.Data[col_name].tolist()]
        for cellLi_index in range(len(cellLi)):
            if cellLi[cellLi_index] == ["."]:
                emptyIndexLi.append(cellLi_index)
                continue
            mergeIndexLi = []  # cellLi_index that needs to be merged
            for k, li in clusterIndexDir.items():
                if set(cellLi[cellLi_index]).intersection(set(li)) != set():
                    mergeIndexLi.append(k)
            if len(mergeIndexLi) == 0:
                clusterIndexDir[str(cellLi_index)] = cellLi[cellLi_index]
            else:
                merge_li = []  # The merged list
                merge_li.extend(cellLi[cellLi_index])
                name = "-".join(mergeIndexLi) + "-" + str(cellLi_index)
                for i in mergeIndexLi:
                    merge_li.extend(clusterIndexDir[i])
                    del clusterIndexDir[i]
                clusterIndexDir[name] = list(set(merge_li))
        return clusterIndexDir, emptyIndexLi

    def mergeData(self, col_name):
        """
        Merge rows with duplicate elements based on column index
        :param col_name: Column Index
        :return: The merged data frame
        """
        clusterIndexDir, emptyIndexLi = self.cluster_li(col_name)
        outDf = pd.DataFrame(columns=self.Data.columns)
        for k, v in clusterIndexDir.items():
            index_li = k.split("-")
            selected_rows = self.Data.iloc[index_li]
            # Pandas merges each row of the data frame, connecting each element with ",". If a cell is ".", the cell is
            # ignored, but if every row cell in a column is ".", the "." is retained.
            merged_row = selected_rows.apply(lambda x: ','.join([str(val) for val in x if val != "."]) if any(val != "." for val in x) else ".", axis=0)
            new_df = pd.DataFrame([merged_row], columns=self.Data.columns)
            new_df = new_df.applymap(lambda x: ','.join(set(x.split(','))))  # Deduplication
            outDf = pd.concat([outDf, new_df], axis=0)
        selected_rows = self.Data.iloc[emptyIndexLi]
        outDf = pd.concat([outDf, selected_rows], axis=0)
        return outDf

    def mergeAllData(self):
        """
        Traverse the column indices of the index and merge rows with duplicate elements
        :return: The merged data frame
        """
        tmp = copy.deepcopy(self.Data)
        col_names = list(self.Data.columns)
        for col_name in col_names:
            self.Data = self.mergeData(col_name)
            # Method for reindexing a DataFrame, replacing the original index with the default integer index (0, 1, 2, ...)
            self.Data = self.Data.reset_index(drop=True)
        outDf = copy.deepcopy(self.Data)
        self.Data = copy.deepcopy(tmp)
        del tmp
        return outDf

    def SaveExcel(self, outExcel):
        self.CNSIndexData.to_excel(outExcel, index=False)
        self.CNSIndexData.to_csv(f"{outExcel}.csv", sep='\t', index=False)

    def Establish_Gene_map(self):
        anchorData = pd.read_csv(self.anchorFile, sep='\t', comment='#', usecols=[0, 1], names=[self.que, self.ref])
        queAllGeneData = pd.read_csv(self.queGeneBed, header=None, sep='\t', names=["chrID", "start", "end", "gene", "score", "strand"])
        refAllGeneData = pd.read_csv(self.refGeneBed, header=None, sep='\t', names=["chrID", "start", "end", "gene", "score", "strand"])
        queAllGeneLi = queAllGeneData["gene"].tolist()
        queMapGeneLi = anchorData[self.que].tolist()
        queUnMapGeneLi = list(set(queAllGeneLi) - set(queMapGeneLi))
        queUnMapGeneData = pd.DataFrame({
            self.que: queUnMapGeneLi,
            self.ref: ["."] * len(queUnMapGeneLi)
        })
        refAllGeneLi = refAllGeneData["gene"].tolist()
        refMapGeneLi = anchorData[self.ref].tolist()
        refUnMapGeneLi = list(set(refAllGeneLi) - set(refMapGeneLi))
        refUnMapGeneData = pd.DataFrame({
            self.que: ["."] * len(refUnMapGeneLi),
            self.ref: refUnMapGeneLi
        })
        self.Data = pd.concat([anchorData, queUnMapGeneData, refUnMapGeneData], ignore_index=True)
        self.Data.to_excel(os.path.join(self.IndexDir, "merge.xlsx"), index=False)
        self.Data.to_csv(os.path.join(self.IndexDir, "merge.csv"), sep='\t', index=False)
        BasePan.pymkdir(os.path.join(self.IndexDir, self.que))
        outExcel = os.path.join(self.IndexDir, self.que, f"{self.que}_{self.ref}.xlsx")
        self.CNSIndexData = self.mergeAllData()
        self.CNSIndexData = self.CNSIndexData[[self.ref, self.que]]  # Swap the positions of the que and ref columns
        self.SaveExcel(outExcel)
        return self.CNSIndexData


class GeneIndexMerge(BasePan):
    """
    1. Import the results of geneIndexAssign
    2. For an index with only one species, find the best match among all species according to diamond. If the hit is in
       the group, merge it; if not, classify the index as Un.
    3. For those > 1, the number of overlaps with other indices is calculated. If overlap >= non overlap, they are
       merged.
    """

    def __init__(self, group_name, group_data, speciesList):
        self.group_name = group_name
        self.group_data = group_data.reset_index(drop=True)
        self.speciesList = speciesList

    @staticmethod
    def merge_cells(cell1, cell2):
        if cell1 == '.' and cell2 == '.':
            return '.'
        elif cell1 == '.':
            return cell2
        elif cell2 == '.':
            return cell1
        else:
            merged_list = list(set(cell1.split(',') + cell2.split(',')))
            return ','.join(sorted(merged_list))

    @staticmethod
    def count_non_dot(row):
        return sum(1 for value in row if value != '.')

    def get_singleSpe(self, row_series):
        for i in self.speciesList:
            if row_series[i] != ".":
                return i, row_series[i]
        else:
            raise Exception("all is .")

    @staticmethod
    def diamond_map_rm(data, target, target_gene, species_li, diamond_dict):
        align_info_li = []
        for i in species_li:
            if i != target:
                map_gene_li = diamond_dict[f"{target}.{i}"].get(target_gene, [])
                map_gene_list = [(x[0], x[1], i) for x in map_gene_li]
                align_info_li.extend(map_gene_list)
        index_li = []
        if len(align_info_li) > 0:
            max_value = max(x[1] for x in align_info_li)
            max_align_info_li = [x for x in align_info_li if x[1] == max_value]
            for align_info in max_align_info_li:
                ref_gene, score, spe = align_info
                contains_Series = data[spe].apply(lambda x: ref_gene in x.split(','))
                if any(contains_Series):
                    contains_index_li = contains_Series[contains_Series].index.tolist()
                    index_li.extend(contains_index_li)
                else:
                    continue
        return list(set(index_li))

    def handle_singleSpe_index(self, diamond_dict):
        """
        It is necessary to ensure that the row index of group_df is unique and the process remains unchanged
        Args:
            diamond_dict: {"A.B": {A_gene1:[(B_gene1, score1), ..., (B_gene2, score2)]}}
        Returns:
            Un_dict: Genes that need to be classified as Un.
            group_df: gene index
        """
        Un_dict = {}
        group_df = copy.deepcopy(self.group_data)
        group_df = group_df.reset_index(drop=True)  # Reset Index
        group_df['num'] = group_df[self.speciesList].apply(self.count_non_dot, axis=1)
        singleSpe_index_li = group_df[group_df['num'] == 1].index.tolist()
        group_df = group_df.drop(columns=['num'])
        for singleSpe_index in singleSpe_index_li:
            rm_row = group_df.loc[singleSpe_index]
            # Maybe because the previous loop added to the row where num=1
            if sum(1 for index, value in rm_row.items() if index not in ["Group", "Index"] and value != ".") == 1:
                group_df = group_df.drop(singleSpe_index, axis=0)
                spe, geneID = self.get_singleSpe(rm_row)
                contains_index_li = self.diamond_map_rm(group_df, spe, geneID, self.speciesList, diamond_dict)
                if len(contains_index_li) == 0:
                    # Put it into group_df and Do not change the index
                    group_df = pd.concat([group_df, rm_row.to_frame().transpose()], axis=0)
                else:
                    for con_index in contains_index_li:
                        group_df.loc[con_index] = group_df.loc[con_index].combine(rm_row, self.merge_cells)
            else:
                pass
        group_df = group_df.reset_index(drop=True)  # Reset Index
        single_spe_df = group_df[(group_df[self.speciesList] != '.').sum(axis=1) == 1]
        group_df = group_df[(group_df[self.speciesList] != '.').sum(axis=1) != 1]
        for _, row in single_spe_df.iterrows():
            spe, geneID = self.get_singleSpe(row)
            Un_dict.setdefault(spe, []).append(geneID)
        return Un_dict, group_df

    def merge_similar_index(self, group_df):
        if len(group_df) != 0:
            rm_data_li = []
            while True:
                group_df = group_df.reset_index(drop=True)
                if group_df.shape[0] == 1:
                    if len(rm_data_li) != 0:
                        tmp_data = pd.concat(rm_data_li, ignore_index=True, axis=0)
                        group_df = pd.concat([group_df, tmp_data], ignore_index=True, axis=0)
                    break
                group_df['dot_count'] = group_df.apply(lambda row: (row == '.').sum(), axis=1)
                group_df = group_df.sort_values(by='dot_count')
                Is_add = False
                max_dot_index = group_df['dot_count'].idxmax()
                group_df = group_df.drop(columns=['dot_count'])
                rm_row = group_df.loc[max_dot_index]
                group_df = group_df.drop(max_dot_index, axis=0)
                for df_index, row in group_df.iterrows():
                    common_count, diff_count = 0, 0
                    for rm_row_col, cell in rm_row.items():
                        if cell == "." or rm_row_col == "Group" or rm_row_col == "Index":
                            continue
                        if set(cell.split(',')) & set(row[rm_row_col].split(',')):
                            common_count += 1
                        else:
                            diff_count += 1
                    if common_count >= diff_count:
                        group_df.loc[df_index] = group_df.loc[df_index].combine(rm_row, self.merge_cells)
                        Is_add = True
                if not Is_add:
                    rm_data_li.append(rm_row.to_frame().transpose())
                else:
                    if len(rm_data_li) != 0:
                        tmp_data = pd.concat(rm_data_li, ignore_index=True, axis=0)
                        group_df = pd.concat([group_df, tmp_data], ignore_index=True, axis=0)
                        rm_data_li = []
        return group_df

    def gene_index_merge(self, diamond_dict):
        Un_dict, group_df = self.handle_singleSpe_index(diamond_dict)
        group_df = self.merge_similar_index(group_df)
        df_li = []
        n = 0
        for sep, gene_li in Un_dict.items():
            n += 1
            df_li.append(pd.DataFrame({
                "Group": [f"{self.group_name}.Un"],
                "Index": [f"{self.group_name}.Un.{n}"],
                sep: [",".join(gene_li)],
            }))
        if len(df_li) > 0:
            tmp_data = pd.concat(df_li, ignore_index=True, axis=0)
            group_df = pd.concat([group_df, tmp_data], ignore_index=True, axis=0).fillna('.')
        return group_df


class ParseTree:
    """
    # Tree Grouping
        1. First, determine whether the tree needs to be split. The basis is that the number of species with paralogous
           genes is greater than the cutoff.
        2. If splitting is required, the tree is split step by step from the root node. After the splitting, the
           subtrees are pruned. That is, there cannot be a single-node subtree. If there is a single-node subtree, the
           single node is assigned to Un.
        3. Determine again whether to continue, if necessary, repeat the above operation. If not, output.
    """
    cutoff = 8  # Determine whether >= cutoff
    smallest_cutoff = 1

    def __init__(self, row, species_li, tree_file):
        self.row = row
        self.species_li = species_li
        self.tree_file = tree_file
        self.gene2spe_dict, self.paralogous_spe_num = self.get_dict()

    @classmethod
    def set_cutoff(cls, cutoff_):
        cls.cutoff = cutoff_

    def get_dict(self):
        """
        This function is used to generate the species corresponding to each gene in each row, and record how many
        species have paralogous
        Returns:
            gene2spe_dict: The dictionary key is gene and the value is specie
            paralogous_spe_num: The number of paralogous species
        """
        paralogous_spe_num = 0
        gene2spe_dict = {}
        for specie in self.species_li:
            if self.row[specie] != ".":
                genes = self.row[specie].split(",")
                if len(genes) > 1:
                    paralogous_spe_num += 1
                for gene in genes:
                    gene2spe_dict[gene] = specie
            else:
                pass
        return gene2spe_dict, paralogous_spe_num

    @staticmethod
    def calculate_distances(tree):
        """
        Calculate the distance between each node and the root node
        Args:
            tree: Tree Object
        Returns:
            Each node has a tree object with the dist_from_root attribute
        """
        out_tree = copy.deepcopy(tree)
        root_node = out_tree.get_tree_root()
        for node in out_tree.traverse():
            distance = tree.get_distance(root_node, node)
            if not hasattr(node, 'dist_from_root'):
                setattr(node, 'dist_from_root', distance)
        return out_tree

    @staticmethod
    def is_tree_empty(tree) -> bool:
        return len(tree.children) == 0

    def parse_tree(self):
        """
        Rename each node and calculate the distance of each node
        Returns:
            tree: new Tree Object
        """
        tree = Tree(self.tree_file)
        name_n = 0
        for node in tree.traverse():
            if node.is_leaf():
                node.name = node.name.split("_longest_pep_")[1]
            else:
                name_n += 1
                node.name = f"node_{name_n}"
        # Create a dist_from_root attribute for each node to record the distance between the node and the root
        tree = self.calculate_distances(tree)
        return tree

    def judge_tree(self, tree) -> bool:
        """
        The decision of whether the tree needs to be split is based on whether the number of species with paralogous
        genes in the tree is greater than the cutoff.
        Args:
            tree: Tree Object
        Returns:
            bool
        """
        if self.is_tree_empty(tree):
            return False
        leaves_nodes = tree.get_leaf_names()
        spe_gene_dict = {}
        for i in leaves_nodes:
            spe_gene_dict.setdefault(self.gene2spe_dict[i], set()).add(i)
        if len([v for v in spe_gene_dict.values() if len(v) > 1]) >= self.cutoff:
            return True
        else:
            return False

    def remove_small_subtree(self, tree):
        """
        Input a tree. If the tree has a subtree with less than or equal to smallest_cutoff leaf nodes, remove the subtree.
        Args:
            tree:

        Returns:

        """
        input_tree = copy.deepcopy(tree)
        rm_nodes_li = []
        while True:
            out_nodes_li = []
            tmp_rm_nodes_li = []  # Store the nodes that need to be rm in a loop
            for child in input_tree.children:
                child_deep = copy.deepcopy(child)
                child_leaves = child_deep.get_leaf_names()
                if len(child_leaves) <= self.smallest_cutoff:
                    rm_nodes_li.extend(child_leaves)
                    tmp_rm_nodes_li.extend(child_leaves)
                else:
                    out_nodes_li.extend(child_leaves)
            if len(out_nodes_li) != 0:
                input_tree.prune(out_nodes_li, preserve_branch_length=True)
                if len(tmp_rm_nodes_li) == 0:
                    break
            else:
                input_tree = Tree()  # Empty tree object
                break
        rm_nodes_set = set(rm_nodes_li)
        return input_tree, rm_nodes_set

    def split_tree(self, tree):
        tree_deep = copy.deepcopy(tree)
        tmp_can_split_tree_li = []
        tmp_un_split_tree_li = []
        if self.judge_tree(tree_deep):
            for child in tree_deep.children:
                child_deep = copy.deepcopy(child)
                if self.judge_tree(child_deep):
                    tmp_can_split_tree_li.append(child_deep)
                else:
                    tmp_un_split_tree_li.append(child_deep)
        else:
            tmp_un_split_tree_li.append(tree_deep)
        return tmp_can_split_tree_li, tmp_un_split_tree_li

    def _tree_to_dataframe(self, tree):
        """
        Output the leaf_node in the tree object into a row of data frame
        Args:
            tree: tree object
        Returns:
            One row of a pd.DataFrame
        """
        dict_ = {}
        dict_data = {}
        leaf_names = tree.get_leaf_names()
        # print(leaf_names)
        for i in leaf_names:
            spe = self.gene2spe_dict[i]
            dict_.setdefault(spe, set()).add(i)
        for i in self.species_li:
            if i not in dict_:
                dict_data[i] = "."
            else:
                dict_data[i] = ",".join(dict_[i])
        return pd.DataFrame([dict_data])

    def run(self):
        if self.paralogous_spe_num >= self.cutoff:  # Further classification is needed using evolutionary trees
            raw_tree = self.parse_tree()
            tree = copy.deepcopy(raw_tree)
            leaf_to_keep_li = list(self.gene2spe_dict.keys())
            tree.prune(leaf_to_keep_li, preserve_branch_length=True)
            rm_nodes_set = set()
            can_split_tree_li = [tree]
            un_split_tree_li = []
            while True:
                if len(can_split_tree_li) != 0:
                    tmp_li = []
                    for i in can_split_tree_li:
                        tmp_out_tree, tmp_rm_nodes_set = self.remove_small_subtree(i)
                        rm_nodes_set = rm_nodes_set.union(tmp_rm_nodes_set)
                        if self.judge_tree(tmp_out_tree):
                            tmp_can_split_tree_li, tmp_un_split_tree_li = self.split_tree(tmp_out_tree)
                            un_split_tree_li.extend(tmp_un_split_tree_li)
                            tmp_li.extend(tmp_can_split_tree_li)
                        else:
                            if not self.is_tree_empty(tmp_out_tree):
                                un_split_tree_li.append(tmp_out_tree)
                    can_split_tree_li, tmp_li = tmp_li, can_split_tree_li
                    del tmp_li
                else:
                    break
            data_li = []
            for un_split_tree in un_split_tree_li:
                tmp_out_tree, tmp_rm_nodes_set = self.remove_small_subtree(un_split_tree)
                rm_nodes_set = rm_nodes_set.union(tmp_rm_nodes_set)
                if not self.is_tree_empty(tmp_out_tree):
                    data_li.append(self._tree_to_dataframe(tmp_out_tree))
            if len(data_li) != 0:
                result = pd.concat(data_li, axis=0, ignore_index=True)
                df = pd.DataFrame({
                    'Group': [self.row["Group"]] * len(result),
                    'Index': [f'{self.row["Index"]}.tree_{i+1}' for i in range(len(result))]
                })
                result = pd.concat([df, result], axis=1)
            else:
                result = pd.DataFrame(columns=list(self.row.index))

            # Processing rm node
            if len(rm_nodes_set) != 0:
                rm_nodes_dict_ = {}
                rm_nodes_dict_data = {'Group': self.row["Group"],
                                      'Index': "{}.tree_Un".format(self.row["Index"])}
                for i in rm_nodes_set:
                    spe = self.gene2spe_dict[i]
                    rm_nodes_dict_.setdefault(spe, set()).add(i)
                for i in self.species_li:
                    if i not in rm_nodes_dict_:
                        rm_nodes_dict_data[i] = "."
                    else:
                        rm_nodes_dict_data[i] = ",".join(rm_nodes_dict_[i])
                rm_nodes_df = pd.DataFrame([rm_nodes_dict_data])
                del rm_nodes_dict_, rm_nodes_dict_data
                # merge
                result = pd.concat([result, rm_nodes_df], axis=0, ignore_index=True)

            return result
        else:
            return pd.DataFrame(self.row).transpose()

