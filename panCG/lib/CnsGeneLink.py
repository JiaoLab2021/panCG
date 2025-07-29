import copy
import pandas as pd
from panCG.lib.base import BasePan


class CnsIndexRead(BasePan):
    recall_cds = "(recall_cds)"
    recall_cns = "(recall_cns)"
    recall_nonCE = "(recall_nonCE)"

    def __init__(self, filter_cns_index):
        self.cns_index_data = pd.read_csv(filter_cns_index, sep="\t")

    def ToBed(self, key_col, value_col):
        col_dict = self.cns_index_data.set_index(key_col)[value_col].to_dict()
        dict_li = []
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
                    chrID, start, end = self.match_CNS_coordinate(CNS)
                    row_dict = {"chrID": chrID, "start": start, "end": end, "cns_index": cns_index}
                    dict_li.append(row_dict)
        return pd.DataFrame(dict_li)

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

    def get_valid_cns_count(self, species_li):
        data = copy.deepcopy(self.cns_index_data)
        data["valid_cns_count"] = data[species_li].apply(self.count_valid_cns, axis=1)
        valid_cns_count_df = data[["Index", "valid_cns_count"]]
        valid_cns_count_df.columns = ["cns_index", "valid_cns_count"]
        return valid_cns_count_df


class GeneIndexRead(BasePan):
    def __init__(self, gene_index_file):
        self.gene_index_data = pd.read_csv(gene_index_file, sep="\t")

    def ToDict(self, key_col, value_col):
        out_dict = {}
        col_dict = self.gene_index_data.set_index(key_col)[value_col].to_dict()
        for gene_index, genes in col_dict.items():
            if genes == ".":
                continue
            gene_li = genes.strip().split(",")
            for gene in gene_li:
                out_dict.setdefault(gene, set()).add(gene_index)
        return out_dict

    def get_gene_count(self, species_li):
        data = copy.deepcopy(self.gene_index_data)
        data["gene_count"] = (data[species_li] != ".").sum(axis=1)
        gene_count_df = data[["Index", "gene_count"]]
        gene_count_df.columns = ["gene_index", "gene_count"]
        return gene_count_df


class LinkCnsGeneIndex:
    def __init__(self, cns_anno_file):
        self.cns_anno_file = cns_anno_file

    def parse_anno(self, anno_dict):
        data = pd.read_csv(self.cns_anno_file, sep="\t")
        data = data.rename(columns={'V4': 'cns_index'})
        data["start"] = data["start"] - 1
        data["gene_index"] = data["transcriptId"].apply(lambda x: ",".join(anno_dict[x]))
        # data = data[data["annotation"] != "Distal_Intergenic"]
        data['annotation'] = data['annotation'].str.replace(r'_\(.*?\)', '', regex=True).str.strip()
        data['annotation'] = data['annotation'].str.replace("'", "", regex=False)
        return data

