import pandas as pd
from scipy.stats import mannwhitneyu


class GenePavAsso:
    min_species = 5

    def __init__(self, pan_gene_file, phenotype_file, out_file):
        self.pan_gene_file = pan_gene_file
        self.phenotype_data = self.parse_phenotype(phenotype_file)
        self.out_file = out_file
    
    @classmethod
    def set_min_species(cls, value_):
        cls.min_species = value_

    @staticmethod
    def parse_phenotype(phenotype_file):
        df = pd.read_csv(phenotype_file, sep="\t", usecols=[0, 1], header=None, names=["Species", "Trait"])
        return df
    
    def gene_pav_asso(self):
        index_df = pd.read_csv(self.pan_gene_file, sep="\t")
        index_df.set_index("Index", inplace=True)
        index_df.drop(["Group"], axis=1, inplace=True)
        out_df_li = []
        for index, row in index_df.iterrows():
            if (row == ".").sum() < self.min_species or (row == ".").sum() > len(row) - self.min_species:
                p_dict = {"Index": index, "p_value": 1.00, "type": "."}
            else:
                row_df = row.to_frame(name='gene_id')
                row_df = row_df.reset_index().rename(columns={'index': 'Species'})
                merge_df = pd.merge(self.phenotype_data, row_df, on="Species", how='left')
                merge_df["x"] = merge_df["gene_id"].apply(lambda x: "with" if x != "." else "without")
                # Perform the two-sided Mann-Whitney U test
                group_A = merge_df[merge_df['x'] == 'with']['Trait']
                group_B = merge_df[merge_df['x'] == 'without']['Trait']
                stat, p = mannwhitneyu(group_A, group_B, alternative='two-sided')
                if group_A.mean() >= group_B.mean():
                    type_ = "with"
                else:
                    type_ = "without"
                p_dict = {"Index": index, "p_value": p, "type": type_}
            out_df_li.append(pd.DataFrame([p_dict]))
        result = pd.concat(out_df_li, axis=0, ignore_index=True)
        result.to_csv(self.out_file, sep='\t', header=True, index=False)
