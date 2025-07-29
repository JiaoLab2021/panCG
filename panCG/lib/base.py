import yaml
import os
import time
import re
import subprocess
from Bio import SeqIO
import pandas as pd
from datetime import datetime
from ete3 import Tree
import shutil
import resource
from functools import wraps

class BasePan:
    @staticmethod
    def read_yaml(yamlFile):
        f = open(yamlFile, "r")
        data = yaml.safe_load(f)
        f.close()
        return data

    @staticmethod
    def write_yaml(data, yamlFile):
        with open(yamlFile, "w") as yaml_file:
            yaml.dump(data, yaml_file, default_flow_style=False, allow_unicode=True)

    @staticmethod
    def mkdirNewDir(workDir, prefix=""):
        now_time = datetime.now().strftime('%Y-%m%d-%H%M%S')
        new_dir = os.path.join(workDir, prefix + now_time)
        while os.path.exists(new_dir):
            time.sleep(1)
            now_time = datetime.now().strftime('%Y-%m%d-%H%M%S')
            new_dir = os.path.join(workDir, prefix + now_time)
        os.makedirs(new_dir)
        return os.path.abspath(new_dir)

    @staticmethod
    def pymkdir(path):
        if not os.path.exists(path):
            os.makedirs(path)  # Creating multi-level directory

    @staticmethod
    def touch_file(file_path):
        if not os.path.exists(file_path):
            with open(file_path, 'w') as f:
                pass

    @staticmethod
    def write_to(file, content):
        f = open(file, "w")
        f.write(content)
        f.close()

    @staticmethod
    def cmd_linux(cmd):
        subprocess.run(cmd, shell=True, executable="/bin/bash", close_fds=True)

    @staticmethod
    def merge_files(output_file, *args, uniq=True):
        """
        This function merges the contents of multiple files into an output file,
        with the same result as the cat command.
        Args:
            output_file: Output after merging
            *args: Files that need to be merged
            uniq: Indicates whether uniq operation needs to be performed on the merged files. [default: True]
        Returns:
            None
        """
        unique_lines = set()
        with open(output_file, 'w') as outfile:
            for fname in args:
                if os.path.isfile(fname):
                    with open(fname, 'r') as infile:
                        for line in infile:
                            if uniq:
                                unique_lines.add(line)
                            else:
                                outfile.write(line)
                else:
                    raise Exception(f"{fname} is not exist")
            if uniq:
                for line in unique_lines:
                    outfile.write(line)

    @staticmethod
    def match_CNS_coordinate(CNSname) -> tuple:
        pattern = r'^(.*?):(\d+)-(\d+)'
        match = re.match(pattern, CNSname)
        if match:
            chrID, start, end = match.groups()
            return chrID, int(start), int(end)
        else:
            raise Exception(
                f"The format of {CNSname} is incorrect. It must be chrID:start-end. ':' cannot appear in chrID")
    
    @staticmethod
    def get_species_li(species_tree, target_species):
        tree = Tree(species_tree)
        target_node = None
        for node in tree.traverse():
            if node.name == target_species:
                target_node = node
                break
        distances = []
        for node in tree.traverse():
            if node != target_node and node.is_leaf():
                distance = target_node.get_distance(node)
                distances.append((node.name, distance))
        sorted_distances = sorted(distances, key=lambda x: x[1])
        species_li = [target_species]
        species_li.extend([i[0] for i in sorted_distances])
        return species_li

    @staticmethod
    def check_cmd(cmd: str) -> str:
        """
        用于检查命令是否存在在当前PATH中
        Args:
            cmd: 命令，如blastp
        Returns:
            path：命令所在的绝对路径
        """
        path_ = shutil.which(cmd)
        if path_:
            return path_
        else:
            raise FileNotFoundError(f"{cmd} not found in PATH")



class TimerDecorator:
    def __init__(self, logger):
        self.logger = logger

    def __call__(self, func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            start_resources = resource.getrusage(resource.RUSAGE_SELF)
            start_time = time.time()
            result = func(*args, **kwargs)
            end_resources = resource.getrusage(resource.RUSAGE_SELF)
            end_time = time.time()
            user_time = end_resources.ru_utime - start_resources.ru_utime
            sys_time = end_resources.ru_stime - start_resources.ru_stime
            self.logger.info(f"{func.__name__} User time: {user_time:.2f} seconds")
            self.logger.info(f"{func.__name__} System time: {sys_time:.2f} seconds")
            self.logger.info(f"{func.__name__} Wall clock time: {end_time - start_time:.2f} seconds")
            return result
        return wrapper


# class TimerDecorator:
#     def __init__(self, logger):
#         self.logger = logger
#
#     def __call__(self, func):
#         def wrapper(*args, **kwargs):
#             start_resources = resource.getrusage(resource.RUSAGE_SELF)
#             start_time = time.time()
#             result = func(*args, **kwargs)
#             end_resources = resource.getrusage(resource.RUSAGE_SELF)
#             end_time = time.time()
#             user_time = end_resources.ru_utime - start_resources.ru_utime
#             sys_time = end_resources.ru_stime - start_resources.ru_stime
#             self.logger.info(f"{func.__name__} User time: {user_time:.2f} seconds")
#             self.logger.info(f"{func.__name__} System time: {sys_time:.2f} seconds")
#             self.logger.info(f"{func.__name__} Wall clock time: {end_time - start_time:.2f} seconds")
#             return result
#         return wrapper


def subtraction_Dict(big_dict, sub_dict):
    """ Calculate the difference between two dictionaries """
    result_dict = {}
    for key in list(big_dict.keys()):
        result_dict[key] = list(set(big_dict[key]) - set(sub_dict.get(key, [])))
    return result_dict


def match_CNS_coordinate(CNSname):
    pattern = r'^(.*?):(\d+)-(\d+)'
    match = re.match(pattern, CNSname)
    if match:
        chrID, start, end = match.groups()
        return chrID, start, end
    else:
        raise Exception(f"The format of {CNSname} is incorrect. It must be chrID:start-end. ':' cannot appear in chrID")


def DetermineReference(SpeCnsDict, speciesList):
    ref = ""
    for i in speciesList:
        if SpeCnsDict[i]:
            ref = i
            break
    if ref != "":
        return ref
    else:
        raise Exception("No reference found")


class GffAnno:
    def __init__(self, logger, gffFile, fastaFile):
        self.logger = logger
        self.gffFile = gffFile
        self.genomeFile = fastaFile
        self.seq_dict = self.parseFastaFile()
        self.mRNA_dict, self.gene_line_dict, self.mRNA_line_dict, self.mRNA_gff_dict = self.parseGffFile()

    def parseFastaFile(self):
        seq_dict = {}
        records = SeqIO.parse(self.genomeFile, "fasta")
        for record in records:
            seq_dict[record.id] = record.seq
        return seq_dict

    def getSeq(self, chrID, start, end):
        """ bed format coordinates """
        if chrID not in self.seq_dict:
            raise ValueError(f"{chrID} not in fasta.")
        seq = self.seq_dict[chrID][start: end]
        return seq

    def get_mRNA(self):
        """
        Returns:
            mRNA_dict: {gene_id: [tran_id, tran_id], ...}
            gene_line_dict: {gene_id: "line"}
            mRNA_line_dict: {tran_id: "line"}
        """
        mRNA_dict = {}
        gene_line_dict = {}
        mRNA_line_dict = {}
        with open(self.gffFile, "r") as f:
            while True:
                line = f.readline()
                if not line:
                    break
                if line.startswith("#") or not line.strip():
                    continue
                line_li = line.strip().split("\t")
                if line_li[2] == "mRNA":
                    gene_id = re.findall(r"Parent=([^;$]+)", line_li[8])[0]
                    mRNA_id = re.findall(r"ID=([^;$]+)", line_li[8])[0]
                    mRNA_dict.setdefault(gene_id, []).append(mRNA_id)
                    mRNA_line_dict[mRNA_id] = line
                    if mRNA_id in mRNA_dict:
                        raise ValueError(f"{mRNA_id} duplicate in gff file")
                elif line_li[2] == "gene":
                    gene_id = re.findall(r"ID=([^;$]+)", line_li[8])[0]
                    if gene_id in gene_line_dict:
                        raise ValueError(f"{gene_id} duplicate in gff file")
                    gene_line_dict[gene_id] = line
        return mRNA_dict, gene_line_dict, mRNA_line_dict

    def parseGffFile(self):
        mRNA_dict, gene_line_dict, mRNA_line_dict = self.get_mRNA()
        mRNA_gff_dict = {mRNA_id: [] for mRNA_li in mRNA_dict.values() for mRNA_id in mRNA_li}
        with open(self.gffFile, "r") as f:
            while True:
                line = f.readline()
                if not line:
                    break
                if line.startswith("#") or not line.strip():
                    continue
                line_li = line.strip().split("\t")
                feature_Parent = re.findall(r"Parent=([^;$]+)", line_li[8])
                if len(feature_Parent) != 0 and feature_Parent[0] in mRNA_gff_dict:
                    mRNA_gff_dict[feature_Parent[0]].append(line)
        return mRNA_dict, gene_line_dict, mRNA_line_dict, mRNA_gff_dict

    def extractCds(self):
        """

        Returns:
            cds_dict: {mRNA_id: cds_seq}
        """
        cds_dict = {}
        for gene_id, mRNA_id_li in self.mRNA_dict.items():
            for mRNA_id in mRNA_id_li:
                cds_dict[mRNA_id] = ""
                data = pd.DataFrame([line.strip().split('\t') for line in self.mRNA_gff_dict[mRNA_id]])
                data.columns = ["chrID", "source", "Type", "start", "end", "score", "strand", "phase", "attributes"]
                cds_data = data[data["Type"] == 'CDS']
                cds_data.loc[:, 'start'] = cds_data['start'].astype(int)
                cds_data.loc[:, 'end'] = cds_data['end'].astype(int)
                if len(cds_data) == 0:
                    self.logger.warning(f"No CDS in mRNA({mRNA_id})")
                    del cds_dict[mRNA_id]
                    continue
                if all(i == "+" for i in cds_data["strand"].tolist()):
                    sorted_cds_data = cds_data.sort_values(by="start")
                    for index, row in sorted_cds_data.iterrows():
                        chrID, start, end = row["chrID"], int(row["start"]) - 1, int(row["end"])
                        cds_dict[mRNA_id] += self.getSeq(chrID, start, end)
                elif all(i == "-" for i in cds_data["strand"].tolist()):
                    sorted_cds_data = cds_data.sort_values(by="start", ascending=False)
                    for index, row in sorted_cds_data.iterrows():
                        chrID, start, end = row["chrID"], int(row["start"]) - 1, int(row["end"])
                        cds_dict[mRNA_id] += self.getSeq(chrID, start, end).reverse_complement()
                else:
                    raise ValueError(f"The CDS positive and negative strands of mRNA({mRNA_id}) exist simultaneously")
        return cds_dict

    def extractPep(self):
        pep_dict = {}
        cds_dict = self.extractCds()
        for mRNA_id, cds in cds_dict.items():
            if len(cds) % 3 != 0:
                self.logger.warning(f"{mRNA_id} not a multiple of three")
            pep_dict[mRNA_id] = str(cds.translate())
        return pep_dict

    def longestGff(self, output_gff_file, longest_type="cds"):
        txt_li = []
        cds_dict = self.extractCds()
        if longest_type == "cds":
            for gene_id, mRNA_id_list in self.mRNA_dict.items():
                tmp_mRNA_dict = {mRNA_id: cds_dict.get(mRNA_id, "") for mRNA_id in mRNA_id_list}
                longest_key = max(tmp_mRNA_dict, key=lambda x: len(tmp_mRNA_dict[x]))
                if len(tmp_mRNA_dict[longest_key]) > 0:
                    txt_li.append(self.gene_line_dict[gene_id])
                    txt_li.append(self.mRNA_line_dict[longest_key])
                    txt_li.append("".join(self.mRNA_gff_dict[longest_key]))
        elif longest_type == "mRNA":
            for gene_id, mRNA_id_list in self.mRNA_dict.items():
                tmp_mRNA_dict = {}
                for mRNA_id in mRNA_id_list:
                    if mRNA_id not in cds_dict:
                        tmp_mRNA_dict[mRNA_id] = 0
                    else:
                        line_li = self.mRNA_line_dict[mRNA_id].strip().split("\t")
                        tmp_mRNA_dict[mRNA_id] = int(line_li[4]) - int(line_li[3])
                longest_key = max(tmp_mRNA_dict, key=tmp_mRNA_dict.get)
                if tmp_mRNA_dict[longest_key] > 0:
                    txt_li.append(self.gene_line_dict[gene_id])
                    txt_li.append(self.mRNA_line_dict[longest_key])
                    txt_li.append("".join(self.mRNA_gff_dict[longest_key]))
        txt = ''.join(txt_li)
        fo = open(output_gff_file, "w")
        fo.write(txt)
        fo.close()


class ParallelScheduler:  
    def __init__(self, logger, total_tasks, completed_tasks, progress_block):  
        self.logger = logger  
        self.total_tasks = total_tasks  
        self.completed_tasks = completed_tasks  
        self.progress_block = progress_block  
  
    def make_call_back(self):  
        def track_progress(_):  
            self.completed_tasks.value += 1  
            if self.completed_tasks.value % self.progress_block == 0:  
                self.logger.info("Completed: {:.2f}%".format((self.completed_tasks.value / self.total_tasks) * 100))  
        return track_progress  
  
    def make_error_callback(self):  
        def track_progress(error):  
            self.logger.error(error)  
        return track_progress  


# class MyBedTool:
#     def __init__(self, bed_data):
#         self.bed_data = bed_data
#         self.rename_bed_data = self.check_bed()
#
#     def check_bed(self):
#         if len(self.bed_data.columns) < 3:
#             raise ValueError("The data frame must contain at least three columns")
#         else:
#             df = self.bed_data.rename(columns={self.bed_data.columns[0]: 'chrID',
#                                                self.bed_data.columns[1]: 'start',
#                                                self.bed_data.columns[2]: 'end'})
#         if df[df.columns[0]].dtype != 'object':
#             raise ValueError("The first column is not of type string")
#         if df[df.columns[1]].dtype != 'int64':
#             raise ValueError("The second column is not of type integer")
#         if df[df.columns[2]].dtype != 'int64':
#             raise ValueError("The third column is not of type integer")
#         return df
#
#     def sort_bed(self):
#         df_sorted = self.rename_bed_data.sort_values(by=['chrID', 'start', 'end'])
#         return df_sorted
#
#     def merge_bed(self, maxGap=0):
#         sort_data = self.sort_bed()
#         sort_data['prev_end'] = sort_data['end'].shift(1)
#         sort_data['prev_chr'] = sort_data['chrID'].shift(1)
#         sort_data['is_overlap'] = (sort_data['start'] <= sort_data['prev_end'] + maxGap) & (
#                     sort_data['chrID'] == sort_data['prev_chr'])
#         # Use the cumsum() function to assign cluster labels to adjacent regions
#         sort_data['cluster'] = (sort_data['is_overlap'] == False).cumsum()
#         sort_data.drop(['prev_end', 'prev_chr', 'is_overlap'], axis=1, inplace=True)
#         grouped = sort_data.groupby(['chrID', 'cluster'], sort=False)
#         merge_bed_dict = {'chrID': [], 'start': [], 'end': []}
#         for tmp_name, tmp_group in grouped:
#             merge_bed_dict['chrID'].append(tmp_name[0])
#             merge_bed_dict['start'].append(tmp_group['start'].min())
#             merge_bed_dict['end'].append(tmp_group['end'].max())
#         merge_bed_data = pd.DataFrame(merge_bed_dict)
#         return merge_bed_data


def split_base_group_data(data, group_by, block_size):
    grouped = data.groupby(group_by, sort=False)
    blockData_li = []
    tmp_data_li = []
    for i, (_, group_data) in enumerate(grouped, start=1):
        if i % block_size == 0:
            tmp_data_li.append(group_data)
            blockData_li.append(pd.concat(tmp_data_li, axis=0))
            tmp_data_li = []
        else:
            tmp_data_li.append(group_data)
    else:
        if len(tmp_data_li) != 0:
            blockData_li.append(pd.concat(tmp_data_li, axis=0))
        del tmp_data_li
    return blockData_li
