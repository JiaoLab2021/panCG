import pandas as pd
import os
import multiprocessing
import random

from panCG.lib.base import BasePan
from panCG.lib.base import ParallelScheduler
from panCG.lib.base import split_base_group_data
from panCG.lib.pancns import read_blast_map


class IndexParse:
    minCNSnum = 2

    def __init__(self, data, speciesList, Reference, logger):
        self.data = data.reset_index(drop=True)
        self.speciesList = speciesList
        self.Reference = Reference
        self.logger = logger

    @classmethod
    def set_minCNSnum(cls, value):
        cls.minCNSnum = value

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
    def blastn_map_rm(data, target, target_cns, species_li, blastn_Dict) -> list:
        """
        Traverse species_li and find the index list of the best matching cns after traversal
        Args:
            data: There is cns data in the species_li column
            target: specie of target cns
            target_cns: target cns
            species_li: species_li
            blastn_Dict: blastn dict
        Returns:
            index list
        """
        align_info_li = []
        for i in species_li:
            if i != target:
                map_cns_li = blastn_Dict[f"{target}.{i}"].get(target_cns, [])
                map_cns_list = [(x[0], x[1], i) for x in map_cns_li]
                align_info_li.extend(map_cns_list)
        sorted_align_info_li = sorted(align_info_li, key=lambda x: x[1], reverse=True)
        for align_info in sorted_align_info_li:
            ref_cns, score, spe = align_info
            contains_Series = data[spe].apply(lambda x: ref_cns in x.split(','))
            if any(contains_Series):
                contains_index_li = contains_Series[contains_Series].index.tolist()
                return contains_index_li
            else:
                continue
        return []

    def mergeSmallerIndex(self, blastn_Dict):
        """
        Merge indexes with smaller number of non-reference columns
        1. Take out the small index
        2. if n=1: Directly merge into other indexes through best blastn score
        3. if 1 < n <= cutoff:
            First, determine whether the CNS of the corresponding column exists in other indexes.
            If it does, merge it directly.
            If it does not, take the index with the best blastn score.
        Args:
            blastn_Dict: {"A.B": {A_cns1:[(B_cns1, score1), ..., (B_cns2, score2)]}}
        Returns:
        """
        if sum(1 for col in self.data.columns if any(self.data[col] != '.')) == 1:
            pass
        else:
            tmp_ref_data_li = []
            while True:
                if self.data.shape[0] == 1:
                    break
                min_count_row = (self.data[self.speciesList]
                                 .apply(lambda row: sum(1 for cell in row if cell != '.'), axis=1)
                                 .idxmin()
                                 )  # Returns the row index of the row with the smallest number of non-"."
                rm_Series = self.data.loc[min_count_row]
                self.data = self.data.drop(min_count_row, axis=0)
                remove_cns_dict = {i: rm_Series[i] for i in self.speciesList if rm_Series[i] != "."}
                if remove_cns_dict.get(self.Reference, ".") != ".":
                    tmp_ref_data_li.append(rm_Series.to_frame().transpose())
                    continue
                if len(remove_cns_dict) <= self.minCNSnum:
                    # For less than minCNSnum, if there is overlap, they are merged;
                    # if there is no overlap, the best blastn alignment is used.
                    if len(remove_cns_dict) == 1:
                        target = list(remove_cns_dict.keys())[0]
                        target_cns = remove_cns_dict[target]
                        contains_index_li = self.blastn_map_rm(self.data, target, target_cns,
                                                               self.speciesList, blastn_Dict)
                        if len(contains_index_li) == 0:
                            tmp_ref_data_li.append(rm_Series.to_frame().transpose())
                        else:
                            for con_index in contains_index_li:
                                self.data.loc[con_index] = self.data.loc[con_index].combine(rm_Series, self.merge_cells)
                    else:
                        tmp_dict = {}
                        for index, row in self.data.iterrows():
                            tmp_dict[index] = 0
                            for column in self.speciesList:
                                row_values = [i.strip() for i in row[column].split(',') if i.strip() != "."]
                                series_values = [i.strip() for i in rm_Series[column].split(',') if i.strip() != "."]
                                if set(row_values) & set(series_values):
                                    tmp_dict[index] += 1
                        max_value = max(tmp_dict.values())
                        if max_value != 0:  # There is row identical to rm_Series
                            for k, v in tmp_dict.items():
                                if v == max_value:
                                    self.data.loc[k] = self.data.loc[k].combine(rm_Series, self.merge_cells)
                        else:  # There is no row identical to rm_Series
                            for target in list(remove_cns_dict.keys()):
                                for target_cns in remove_cns_dict[target].split(","):
                                    contains_index_li = self.blastn_map_rm(self.data, target, target_cns,
                                                                           self.speciesList, blastn_Dict)
                                    if len(contains_index_li) > 0:
                                        for con_index in contains_index_li:
                                            self.data.loc[con_index] = (self.data.loc[con_index]
                                                                        .combine(rm_Series, self.merge_cells))
                                        break
                                else:
                                    continue
                                break
                            else:
                                # self.logger.warning("There is no blast hit in this index,
                                # so it is added directly to the index. {}".format(rm_Series.to_dict()))
                                tmp_ref_data_li.append(rm_Series.to_frame().transpose())
                else:
                    self.data = pd.concat([self.data, rm_Series.to_frame().transpose()],
                                          ignore_index=True, axis=0)  # Add back the deleted line
                    break
            if len(tmp_ref_data_li) != 0:
                tmp_ref_data = pd.concat(tmp_ref_data_li, ignore_index=True, axis=0)
                self.data = pd.concat([self.data, tmp_ref_data], ignore_index=True, axis=0)

    def mergeSimilarIndex(self):
        rm_data_li = []
        while True:
            self.data = self.data.reset_index(drop=True)
            if self.data.shape[0] == 1:
                if len(rm_data_li) != 0:
                    tmp_data = pd.concat(rm_data_li, ignore_index=True, axis=0)
                    self.data = pd.concat([self.data, tmp_data], ignore_index=True, axis=0)
                break
            self.data['dot_count'] = self.data.apply(lambda row: (row == '.').sum(), axis=1)
            self.data = self.data.sort_values(by='dot_count')
            Is_add = False
            max_dot_index = self.data['dot_count'].idxmax()
            self.data = self.data.drop(columns=['dot_count'])
            rm_row = self.data.loc[max_dot_index]
            self.data = self.data.drop(max_dot_index, axis=0)
            for df_index, row in self.data.iterrows():
                common_count, diff_count = 0, 0
                for rm_row_col, cell in rm_row.items():
                    if cell == "." or rm_row_col == "Group" or rm_row_col == "Index":
                        continue
                    if set(cell.split(',')) & set(row[rm_row_col].split(',')):
                        common_count += 1
                    else:
                        diff_count += 1
                if common_count >= diff_count:
                    self.data.loc[df_index] = self.data.loc[df_index].combine(rm_row, self.merge_cells)
                    Is_add = True
            if not Is_add:
                rm_data_li.append(rm_row.to_frame().transpose())
            else:
                if len(rm_data_li) != 0:
                    tmp_data = pd.concat(rm_data_li, ignore_index=True, axis=0)
                    self.data = pd.concat([self.data, tmp_data], ignore_index=True, axis=0)
                    rm_data_li = []
        del rm_data_li


def multi_IndexMerge(data, speciesList, Reference, logger, blastn_Dict):
    Grouped = data.groupby('Group', sort=False)
    group_data_li = []
    for _, group_data in Grouped:
        if group_data.shape[0] == 1:
            group_data_li.append(group_data)
        else:
            IndexParser = IndexParse(group_data, speciesList, Reference, logger)
            IndexParser.mergeSmallerIndex(blastn_Dict)
            IndexParser.mergeSimilarIndex()
            group_data_li.append(IndexParser.data)
    return pd.concat(group_data_li, axis=0)


def run_pancnsMerge(logger, config, workDir, Reference, args_dict):
    logger.info("---------------------------------- step 4. Start cnsIndexMerge ... ----------------------------------")
    RefIndexDir = os.path.join(workDir, "Ref_{}_IndexDir".format(Reference))
    cnsIndexAssignFile = os.path.join(RefIndexDir, "Ref.{}.cnsIndexAssign.csv".format(Reference))
    configData = BasePan.read_yaml(config)
    IndexParse.set_minCNSnum(args_dict["min_cns_num"])
    species_tree = configData["species_tree"]
    speciesList = BasePan.get_species_li(species_tree, Reference)
    blastnDir = os.path.join(workDir, "blastnDir")

    # Store all pair-specie filter blastn CNS in the dictionary
    logger.info("Store all pair-specie filter blastn CNS in the dictionary ...")
    blastn_Dict = {}
    for que in speciesList:
        for ref in speciesList:
            if que == ref:
                continue
            blastn_File = os.path.join(blastnDir, f"{que}.{ref}",
                                       f"{que}.{ref}.blastn.halLiftoverFilter.anchors.OriginalFmt.txt")
            blastn_Dict[f"{que}.{ref}"] = read_blast_map(blastn_File)
    logger.info("Finish data loading")

    data = pd.read_table(cnsIndexAssignFile, sep="\t")
    block_size = args_dict["merge_chunk_size"]
    blockData_li = split_base_group_data(data, "Group", block_size)

    random_seed = args_dict["merge_random_seed"]
    random.seed(random_seed)
    # Disrupt the order of the blockData_li list to prevent unbalanced input data during parallel processing
    random.shuffle(blockData_li)

    total_tasks = len(blockData_li)
    threads = args_dict["merge_threads"]
    logger.info("total_tasks: {}".format(total_tasks))
    logger.info("threads: {}".format(threads))
    completed_tasks = multiprocessing.Manager().Value('i', 0)
    progress_block = threads
    parallelScheduler = ParallelScheduler(logger, total_tasks, completed_tasks, progress_block)
    pool = multiprocessing.Pool(processes=threads)
    results = []
    for blockData in blockData_li:
        track_progress = parallelScheduler.make_call_back()
        error_callback = parallelScheduler.make_error_callback()
        result = pool.apply_async(multi_IndexMerge,
                                  args=(blockData, speciesList, Reference, logger, blastn_Dict,),
                                  callback=track_progress, error_callback=error_callback)
        results.append(result)
    pool.close()
    pool.join()

    cnsIndexMergeData = pd.concat([i.get() for i in results], axis=0)
    RefIndexDir = os.path.join(workDir, "Ref_{}_IndexDir".format(Reference))
    cnsIndexMergeFile = os.path.join(RefIndexDir, "Ref.{}.cnsIndexMerge.csv".format(Reference))
    cnsIndexMergeData.to_csv(cnsIndexMergeFile, sep='\t', index=False)

    logger.info("End cnsIndexMerge !!!")
