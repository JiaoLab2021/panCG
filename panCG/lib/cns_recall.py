import pandas as pd
import pyranges as pr
from panCG.lib.base import BasePan


class HalLiftoverMapCNSs(BasePan):
    AverageBwScoreThreshold = None
    effective_rate = 0.7  # The minimum value of effective_len / max(Candidate_ce, exist_cns) during recall
    cns_rate = 0.5  # The minimum value of overlap(Candidate_ce, exist_cns) / min(Candidate_ce, exist_cns),
    # if it exceeds, use exist cns instead

    def __init__(self, group_data, all_halLiftover_dict, speciesList, workDir):
        self.group_data = group_data
        self.all_halLiftover_dict = all_halLiftover_dict
        self.speciesList = speciesList
        self.workDir = workDir

    @classmethod
    def set_AverageBwScoreThreshold(cls, AverageBwScoreThreshold):
        cls.AverageBwScoreThreshold = AverageBwScoreThreshold

    @classmethod
    def set_effective_rate(cls, effective_rate):
        cls.effective_rate = effective_rate

    @classmethod
    def set_cns_rate(cls, cns_rate):
        cls.cns_rate = cns_rate

    def col2bed(self, col_name) -> list:
        cns_li = []
        for cell in self.group_data[col_name].tolist():
            if cell != ".":
                for cns in cell.split(","):
                    chrID, start, end = self.match_CNS_coordinate(cns)
                    cns_li.append([chrID, start, end])
        return cns_li

    @staticmethod
    def get_overlap_len(region_1: list, region_2: list) -> int:
        """
        Calculate the overlap length of two genomic coordinates
        Args:
            region_1: [chr, start, end]
            region_2: [chr, start, end]
        Returns:
            overlap length
        """
        if not region_1[1] <= region_1[2] or not region_2[1] <= region_2[2]:
            raise ValueError("end must be greater than start")
        if region_1[0] == region_2[0]:
            if region_1[1] <= region_2[1]:
                a, b = region_1, region_2
            else:
                a, b = region_2, region_1
            if a[2] <= b[1]:
                return 0
            elif b[1] < a[2] <= b[2]:
                return a[2] - b[1]
            else:
                return b[2] - b[1]
        else:
            return 0

    def refCNS2queCE(self, ref_start, ref_end, que, Candidate_ce_li):
        """

        Args:
            ref_start: in reference
            ref_end: in reference
            que: que
            Candidate_ce_li: in que as [[(chrID, start, end), (score, len)], [(chrID, start, end), (score, len)]]

        Returns:

        """
        CEs_Dict = {"CNS": set(), "CEs": set(), "noCEs": set()}
        data_cns_li = self.col2bed(que)
        for Candidate_ce in Candidate_ce_li:
            Candidate_chrID, Candidate_start, Candidate_end = Candidate_ce[0]
            averageBwScore, effective_len = Candidate_ce[1]
            if effective_len / max(Candidate_end - Candidate_start, ref_end - ref_start) >= self.effective_rate:
                if averageBwScore >= self.AverageBwScoreThreshold:
                    # Determine whether the obtained Candidate_ce needs to be merged with the existing index
                    if len(data_cns_li) > 0:  # This column has cns in other indexes
                        overlap_len_dict = {}
                        for data_cns in data_cns_li:  # Determine whether there is a large overlap with the identified CNS
                            data_cns_name = f"{data_cns[0]}:{data_cns[1]}-{data_cns[2]}"
                            overlap_len_dict[data_cns_name] = self.get_overlap_len(Candidate_ce[0], data_cns)
                        max_overlap_len_key = max(overlap_len_dict, key=overlap_len_dict.get)
                        max_overlap_len = overlap_len_dict[max_overlap_len_key]
                        if max_overlap_len != 0:  # overlap
                            data_cns_chrID, data_cns_start, data_cns_end = self.match_CNS_coordinate(max_overlap_len_key)
                            if max_overlap_len / min(Candidate_end - Candidate_start,
                                                     data_cns_end - data_cns_start) >= self.cns_rate:
                                CEs_Dict["CNS"].add(max_overlap_len_key)
                            else:
                                CEs_Dict["CEs"].add(f"{Candidate_chrID}:{Candidate_start}-{Candidate_end}(CEs)")
                        else:  # not overlap
                            CEs_Dict["CEs"].add(f"{Candidate_chrID}:{Candidate_start}-{Candidate_end}(CEs)")
                    else:  # This column has no cns in other indexes
                        CEs_Dict["CEs"].add(f"{Candidate_chrID}:{Candidate_start}-{Candidate_end}(CEs)")
                else:  # Conservative score does not exceed the threshold
                    CEs_Dict["noCEs"].add(f"{Candidate_chrID}:{Candidate_start}-{Candidate_end}(recall_nonCE)")
            else:
                pass
        ce_set = set()
        ce_set = ce_set.union(CEs_Dict["CNS"], CEs_Dict["CEs"], CEs_Dict["noCEs"])
        return list(ce_set)

    def run(self):
        data_li = []
        for index, row in self.group_data.iterrows():
            reference = ""
            row_dict = row.to_dict()
            Index_Dict = {}  # Store the CNS of each Index after recall
            # Determine the reference genome. By default, the first column is the reference genome.
            for specie in self.speciesList:
                if row_dict[specie] != ".":
                    reference = specie
                    break
            else:
                raise ValueError(f"The line does not exist cns. \n{row_dict}")

            ref_cns_li = row_dict[reference].split(",")
            if len(ref_cns_li) == 1:
                chrID, start, end = self.match_CNS_coordinate(row_dict[reference].split(",")[0])
            else:  # A reference may have multiple cns
                li = [(x[0], x[1], x[2]) for x in (self.match_CNS_coordinate(i) for i in ref_cns_li)]
                li_sorted = sorted(li, key=lambda x: x[2] - x[1], reverse=True)
                chrID, start, end = li_sorted[0]
            for specie, CNSsStr in row_dict.items():
                if specie == "Group" or specie == "Index":  # Skip the first two columns
                    Index_Dict[specie] = CNSsStr
                    continue
                if CNSsStr == ".":
                    # all_halLiftover_dict as [[(chrID, start, end), (score, len)], [(chrID, start, end), (score, len)]]
                    Candidate_ce_li = (self.all_halLiftover_dict[f"{reference}.{specie}"]
                                       .get(f"{chrID}:{start}-{end}", []))
                    if len(Candidate_ce_li) == 0:
                        Index_Dict[specie] = CNSsStr
                    else:
                        ce_li = self.refCNS2queCE(start, end, specie, Candidate_ce_li)
                        if len(ce_li) != 0:
                            Index_Dict[specie] = ",".join(ce_li)
                        else:
                            Index_Dict[specie] = "."
                else:
                    Index_Dict[specie] = CNSsStr
            data_li.append(pd.DataFrame([Index_Dict]))
        return pd.concat(data_li, ignore_index=True, axis=0)


class ReCnsIndexMerge(BasePan):
    def __init__(self, group_data, speciesList):
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

    def indexMerge(self):
        """
        1. Extract the rows with less cns in group_data and delete them in group_data
        2. Traverse other rows to see if they can be merged with other rows (the requirement for merging is that the
           number of columns with the same cns is greater than the number of columns without the same cns)
        3. If there is no merge, store the deleted row in no_merge_li; if there is a merge, merge the data of
           no_merge_li and continue from step 1
        4. The loop stops until no more merges are possible.
        """
        no_merge_li = []
        while True:
            self.group_data = self.group_data.reset_index(drop=True)
            if self.group_data.shape[0] == 1:
                if len(no_merge_li) != 0:
                    no_merge_data = pd.concat(no_merge_li, ignore_index=True, axis=0)
                    self.group_data = pd.concat([self.group_data, no_merge_data], ignore_index=True, axis=0)
                break
            self.group_data['dot_count'] = self.group_data.apply(lambda row: (row == '.').sum(), axis=1)
            self.group_data = self.group_data.sort_values(by='dot_count')
            Is_add = False
            max_dot_index = self.group_data['dot_count'].idxmax()
            self.group_data = self.group_data.drop(columns=['dot_count'])
            rm_row = self.group_data.loc[max_dot_index]
            self.group_data = self.group_data.drop(max_dot_index, axis=0)
            for df_index, row in self.group_data.iterrows():
                common_count, diff_count = 0, 0
                for rm_row_col, cell in rm_row.items():
                    if cell == "." or rm_row_col == "Group" or rm_row_col == "Index":
                        continue
                    if set(cell.split(',')) & set(row[rm_row_col].split(',')):
                        common_count += 1
                    else:
                        diff_count += 1
                if common_count >= diff_count:
                    self.group_data.loc[df_index] = self.group_data.loc[df_index].combine(rm_row, self.merge_cells)
                    Is_add = True
            if not Is_add:
                no_merge_li.append(rm_row.to_frame().transpose())
            else:
                if len(no_merge_li) != 0:
                    tmp_data = pd.concat(no_merge_li, ignore_index=True, axis=0)
                    self.group_data = pd.concat([self.group_data, tmp_data], ignore_index=True, axis=0)
                    no_merge_li = []
        del no_merge_li


class CeCds(BasePan):
    """
    This class handles the overlap of CEs and cds in each column
    """
    def __init__(self, data, col_name, cds_file, out_file):
        self.data = data
        self.col_name = col_name
        self.cds_file = cds_file
        self.out_file = out_file

    def get_ces(self):
        line_li = []
        # txt = ""
        cell_li = self.data[self.col_name].tolist()
        for cell in cell_li:
            if cell == ".":
                continue
            for ce in cell.split(","):
                if "(CEs)" in ce:
                    chrID, start, end = self.match_CNS_coordinate(ce.replace("(CEs)", ""))
                    line_li.append((chrID, start, end, f"{chrID}:{start}-{end}"))
                    # txt += f"{chrID}\t{start}\t{end}\t{chrID}:{start}-{end}\n"
        df = pd.DataFrame(line_li, columns=["Chromosome", "Start", "End", "name"])
        pr_df = pr.PyRanges(df)
        sort_pr_df = pr_df.sort()
        cds_df = pd.read_csv(self.cds_file, sep="\t", header=None, names=["Chromosome", "Start", "End"])
        cds_pr = pr.PyRanges(cds_df)
        intersect_pr = sort_pr_df.join(cds_pr, report_overlap=True)
        intersect_df = intersect_pr.df
        intersect_df = intersect_df[["Chromosome", "Start", "End", "name"]]
        intersect_df.to_csv(self.out_file, sep="\t", header=False, index=False)


class TripleCnsIndexMerge(ReCnsIndexMerge):
    def __init__(self, group_data, speciesList):
        super().__init__(group_data, speciesList)

    @staticmethod
    def merge_region(region_li):
        merged = []
        for region in region_li:
            if not merged or merged[-1][2] < region[1]:
                merged.append(region)
            else:
                merged[-1][2] = max(merged[-1][2], region[2])
        return merged

    @classmethod
    def merge_overlap(cls, cell):
        if cell == ".":
            return cell
        else:
            recall_cns_li, recall_cds_li, cns_li = [], [], []
            values_li = cell.split(",")
            for i in values_li:
                if "(recall_cns)" in i:
                    chrID, start, end = cls.match_CNS_coordinate(i.replace("(recall_cns)", ""))
                    recall_cns_li.append([chrID, start, end])
                elif "(recall_cds)" in i:
                    chrID, start, end = cls.match_CNS_coordinate(i.replace("(recall_cds)", ""))
                    recall_cds_li.append([chrID, start, end])
                else:
                    cns_li.append(i)
            if len(recall_cns_li) >= 2:
                recall_cns_li = sorted(recall_cns_li, key=lambda x: x[1])
                recall_cns_li = cls.merge_region(recall_cns_li)
            if len(recall_cds_li) >= 2:
                recall_cds_li = sorted(recall_cds_li, key=lambda x: x[1])
                recall_cds_li = cls.merge_region(recall_cds_li)
            recall_cns_li = [f"{i[0]}:{i[1]}-{i[2]}(recall_cns)" for i in recall_cns_li]
            recall_cds_li = [f"{i[0]}:{i[1]}-{i[2]}(recall_cds)" for i in recall_cds_li]
            ce_li = [*cns_li, *recall_cns_li, *recall_cds_li]
            return ",".join(ce_li)

    @staticmethod
    def check_merge(series1, series2):
        """
        Determine whether two series need to be merged.
        The merge condition is that in the cells with the same index, the number of recall_cns or recall_cds with
        overlap is dominant.
        Args:
            series1:
            series2:
        Returns:
            bool
        """
        if not series1.index.equals(series2.index):
            raise Exception("The indexes of the two series are inconsistent")
        common_count, diff_count = 0, 0
        for (index1, value1), (index2, value2) in zip(series1.items(), series2.items()):
            if index1 == "Group" or index1 == "Index" or value1 == "." or value2 == ".":
                continue
            value1_li, value2_li = value1.split(","), value2.split(",")
            recall_cns_li_1, recall_cds_li_1, cns_li_1 = [], [], []
            recall_cns_li_2, recall_cds_li_2, cns_li_2 = [], [], []
            for i in value1_li:
                if "(recall_cns)" in i:
                    chrID, start, end = BasePan.match_CNS_coordinate(i.replace("(recall_cns)", ""))
                    recall_cns_li_1.append([chrID, start, end])
                elif "(recall_cds)" in i:
                    chrID, start, end = BasePan.match_CNS_coordinate(i.replace("(recall_cds)", ""))
                    recall_cds_li_1.append([chrID, start, end])
                else:
                    cns_li_1.append(i)
            for i in value2_li:
                if "(recall_cns)" in i:
                    chrID, start, end = BasePan.match_CNS_coordinate(i.replace("(recall_cns)", ""))
                    recall_cns_li_2.append([chrID, start, end])
                elif "(recall_cds)" in i:
                    chrID, start, end = BasePan.match_CNS_coordinate(i.replace("(recall_cds)", ""))
                    recall_cds_li_2.append([chrID, start, end])
                else:
                    cns_li_2.append(i)
            # Determine whether there is overlap between recall_cns_li_1 and recall_cns_li_2 and
            # between recall_cds_li_1 and recall_cds_li_2
            recall_cns_overlap, recall_cds_overlap, cns_overlap = False, False, False
            for i in recall_cns_li_1:
                for j in recall_cns_li_2:
                    overlap_len = HalLiftoverMapCNSs.get_overlap_len(i, j)
                    if overlap_len > 0:
                        recall_cns_overlap = True
                        break
                else:
                    continue
                break
            for i in recall_cds_li_1:
                for j in recall_cds_li_2:
                    overlap_len = HalLiftoverMapCNSs.get_overlap_len(i, j)
                    if overlap_len > 0:
                        recall_cds_overlap = True
                        break
                else:
                    continue
                break
            if set(cns_li_1) & set(cns_li_2):
                cns_overlap = True
            if recall_cns_overlap or recall_cds_overlap or cns_overlap:
                common_count += 1
            else:
                diff_count += 1
        if common_count >= diff_count:
            return True
        else:
            return False

    def indexMerge(self):
        no_merge_li = []
        while True:
            self.group_data = self.group_data.reset_index(drop=True)
            if self.group_data.shape[0] == 1:
                if len(no_merge_li) != 0:
                    no_merge_data = pd.concat(no_merge_li, ignore_index=True, axis=0)
                    self.group_data = pd.concat([self.group_data, no_merge_data], ignore_index=True, axis=0)
                break
            self.group_data['dot_count'] = self.group_data.apply(lambda row: (row == '.').sum(), axis=1)
            self.group_data = self.group_data.sort_values(by='dot_count')
            Is_add = False
            max_dot_index = self.group_data['dot_count'].idxmax()
            self.group_data = self.group_data.drop(columns=['dot_count'])
            rm_row = self.group_data.loc[max_dot_index]
            self.group_data = self.group_data.drop(max_dot_index, axis=0)
            for df_index, row in self.group_data.iterrows():
                if self.check_merge(row, rm_row):
                    self.group_data.loc[df_index] = self.group_data.loc[df_index].combine(rm_row, self.merge_cells)
                    Is_add = True
            if not Is_add:
                no_merge_li.append(rm_row.to_frame().transpose())
            else:
                if len(no_merge_li) != 0:
                    tmp_data = pd.concat(no_merge_li, ignore_index=True, axis=0)
                    self.group_data = pd.concat([self.group_data, tmp_data], ignore_index=True, axis=0)
                    no_merge_li = []
        del no_merge_li

