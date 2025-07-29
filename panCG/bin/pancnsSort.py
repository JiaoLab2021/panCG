import os
import pandas as pd
from panCG.lib.base import BasePan


class Data2SortedBed:
    def __init__(self, data, col_name):
        self.data = data
        self.col_name = col_name

    @staticmethod
    def sort_key(item):
        chr_id, rest = item.split(":")
        start = int(rest.split("-")[0])
        return chr_id, start

    def choose_sort_cns(self, cell):
        cns_li = cell.split(",")
        if len(cns_li) == 1:
            return cell
        else:
            return sorted(cns_li, key=self.sort_key)[0]

    def data2SortedBed(self):
        col_data = (
            self.data
            .loc[:, [self.col_name]]
            .query(f"{self.col_name} != '.'")
        )
        col_data = col_data.assign(**{
            self.col_name: col_data[self.col_name]
                                   .str.replace("(recall_cns)", '', regex=False)
                                   .str.replace("(recall_cds)", '', regex=False)
                                   .str.replace("(recall_nonCE)", '', regex=False)})
        col_data[self.col_name] = col_data[self.col_name].apply(lambda x: self.choose_sort_cns(x))
        result = (
            col_data
            .assign(chr=lambda x: x[self.col_name].str.split(':').str[0],
                    end=lambda x: x[self.col_name].str.split(':').str[1].str.split('-').str[1].astype(int))
            .drop(columns=[self.col_name])
            .reindex(columns=['chr', 'start', 'end'])
            .sort_values(by=['chr', 'start', 'end'])
        )
        return result


def run_pancnsSort(logger, config, workDir, Reference):
    logger.info("----------------------------------- step 6. Start sort PanCNS ... -----------------------------------")
    configData = BasePan.read_yaml(config)
    species_tree = configData["species_tree"]
    speciesList = BasePan.get_species_li(species_tree, Reference)

    RefIndexDir = os.path.join(workDir, "Ref_{}_IndexDir".format(Reference))
    cnsRecall_file = os.path.join(RefIndexDir, "Ref.{}.recall.csv".format(Reference))
    data = pd.read_table(cnsRecall_file)
    data['tmp_index'] = ['tmp{:0{}d}'.format(i, len(str(len(data)))) for i in range(1, len(data) + 1)]
    data = data.set_index('tmp_index')
    
    result_data = pd.DataFrame(columns=data.columns)

    for col in speciesList:  # Group Index
        Data2SortedBeder = Data2SortedBed(data, col)
        df = Data2SortedBeder.data2SortedBed()
        result_data = pd.concat([result_data, data.loc[df.index]], axis=0)
        li = list(set(data.index) - set(df.index))
        data = data.loc[li]
    if len(data) != 0:
        print(data)
        raise Exception("There are rows in the data frame that are not involved in sorting")

    cnsIndexSort_file = os.path.join(RefIndexDir, "Ref.{}.sort.csv".format(Reference))
    result_data.to_csv(cnsIndexSort_file, sep='\t', index=False)

    logger.info("Finish sort PanCNS ...")
    logger.info("The final result is output in {}".format(os.path.realpath(cnsIndexSort_file)))

