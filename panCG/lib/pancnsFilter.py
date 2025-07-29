import os
import re
import pandas as pd


class PanCnsFilter:
    minRow = 3
    recall_cds = "(recall_cds)"
    recall_cns = "(recall_cns)"
    recall_nonCE = "(recall_nonCE)"

    # def __init__(self, PanCns_file):
    #     self.PanCns_file = PanCns_file

    @staticmethod
    def pymkdir(path):
        if not os.path.exists(path):
            os.makedirs(path)

    @staticmethod
    def match_CNS_coordinate(CNSname):
        pattern = r'^(.*?):(\d+)-(\d+)'
        match = re.match(pattern, CNSname)
        if match:
            chrID, start, end = match.groups()
            return chrID, start, end
        else:
            raise Exception(
                f"The format of {CNSname} is incorrect. It must be chrID:start-end. ':' cannot appear in chrID")

    def ToBed(self, data, key_col, value_col):
        col_dict = data.set_index(key_col)[value_col].to_dict()
        tuple_li = []
        for cns_index, cnss in col_dict.items():
            if cnss == ".":
                continue
            cns_li = cnss.strip().split(",")
            for cns in cns_li:
                if self.recall_cds not in cns and self.recall_nonCE not in cns:
                    if self.recall_cns in cns:
                        CNS = cns.replace(self.recall_cns, "")
                    else:
                        CNS = cns
                    tuple_li.append((CNS, cns_index))
        df = pd.DataFrame(tuple_li, columns=['CNS', 'cns_index'])
        df[['chrID', 'start', 'end']] = df['CNS'].apply(self.match_CNS_coordinate).apply(pd.Series)
        df = df[['chrID', 'start', 'end', 'CNS']]
        df = df.drop_duplicates()
        return df

    @staticmethod
    def write_to(file, content):
        f = open(file, "w")
        f.write(content)
        f.close()

    def count_full_cds(self, row):
        count = 0
        for col in row:
            if col != ".":
                elements = col.split(',')
                if all(self.recall_cds in element for element in elements):
                    count += 1
        return count

    def count_valid_cns(self, row):
        """
        Calculate the number of valid CNS, that is, the number of non-cds and non-nonCE
        Args:
            row:
        Returns:
        """
        count = 0
        for col in row:
            if col != ".":
                elements = col.split(',')
                if any((self.recall_cds not in element and self.recall_nonCE not in element) for element in elements):
                    count += 1
        return count

    def count_valid_ce(self, row):
        """
        Calculate the number of valid CEs, that is, the number of non-CEs
        Args:
            row:
        Returns:
        """
        count = 0
        for col in row:
            if col != ".":
                elements = col.split(',')
                if any(self.recall_nonCE not in element for element in elements):
                    count += 1
        return count
