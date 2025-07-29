from Bio import SeqIO
import pandas as pd
import os
import re
import copy
import uuid
import pyranges as pr
from panCG.lib.base import BasePan


def combine_columns(data, *args):
    cols = [data[i] for i in list(args)]
    combined_values_li = [col for col in cols if col != "."] if any(col != "." for col in cols) else []
    combined_values_li = list(set(combined_values_li))
    if len(combined_values_li) != 0:
        return ','.join(combined_values_li)
    else:
        return "."


def read_blast_map(blastFile):
    """
    Read blastn result file and output map dictionary
    Args:
        blastFile: blast result file
    Returns:
        a dictionary type {A_cns1:[(B_cns1, score1), ..., (B_cns2, score2)]}
    """
    Dict = {}
    with open(blastFile, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            line_li = line.strip().split("\t")
            key, value, score = line_li[0], line_li[1], float(line_li[-1])
            Dict.setdefault(key, []).append((value, score))
    return Dict


def read_map(filter_File):
    """
    The input file has only two columns, such as:
    chr1_RagTag:53840-53886 chr1:36854714-36854752
    chr1_RagTag:54015-54132 chr1:36854454-36854501
    chr1_RagTag:54151-54279 chr1:36854320-36854446
    """
    Dict = {}
    with open(filter_File, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            line_li = line.strip().split("\t")
            key, value = line_li[0], line_li[1]
            Dict.setdefault(key, []).append(value)
    return Dict


class CNSsAnchor(BasePan):

    """
    Integrate the results of halLiftover, JCVI and blastn to generate the CNS corresponding
    relationship between que and ref queCNS-refCNS
    """
    halLiftover_cmd = BasePan.check_cmd("halLiftover")
    makeblastdb_cmd = BasePan.check_cmd("makeblastdb")
    blastn_cmd = BasePan.check_cmd("blastn")
    blastn_threads = 6
    jcvi_threads = 6
    overlap_rate = 0.70  # The overlap between the halLiftover area and the CNS / min(halLiftover, CNS) is used for the
    # halLiftover of  the CNS between the two species. If it is greater than this threshold, it is considered that the
    # two CNSs have  passed the halLiftover.
    blastn_evalue = 0.01  # blastn -task blastn-short
    blastn_overlap_rate = 0.5  # To filter the blastn results, add a condition:
    # length of the two CNS alignments / min(CNS1, CNS2) > blastn_overlap_rate
    max_gap = 10

    def __init__(self, HalFile, que, ref, queCNSsBedFile, refCNSsBedFile, queGenomeFile, refGenomeFile, workDir):
        """
        This class integrates the results of halLiftover, jcvi and blastn to realize the CNS map between ref and que.
        Args:
            HalFile: Alignment results of Cactus
            que: query
            ref: reference
            queCNSsBedFile:  The CNS bed file for query must have 4 columns: chrID, start, end, cns_name (chrID:start-end)
            refCNSsBedFile: The CNS bed file for reference must have 4 columns: chrID, start, end, cns_name (chrID:start-end)
            queGenomeFile: query genome file
            refGenomeFile: reference genome file
            workDir: Basic working path
        """
        self.HalFile = HalFile
        self.que = que
        self.ref = ref
        self.queCNSsBedFile = queCNSsBedFile
        self.refCNSsBedFile = refCNSsBedFile
        self.queGenomeFile = queGenomeFile
        self.refGenomeFile = refGenomeFile
        self.workDir = workDir

    @ classmethod
    def modify_blastn_evalue(cls, value):
        cls.blastn_evalue = value

    @ classmethod
    def modify_max_gap(cls, value):
        cls.max_gap = value

    @ classmethod
    def modify_threads(cls, threads_1, threads_2):
        cls.blastn_threads = threads_1
        cls.jcvi_threads = threads_2

    @ classmethod
    def modify_overlap_rate(cls, value):
        cls.overlap_rate = value

    def Run_halLiftover(self):
        # 1. parallel HalLiftover
        halLiftover_bed_file = f"{self.que}.{self.ref}.bed"
        cmd = (f"{self.halLiftover_cmd} {self.HalFile} {self.que} {self.queCNSsBedFile} "
               f"{self.ref} {halLiftover_bed_file}")
        self.cmd_linux(cmd)
        # 2. merge close regions in HalLiftover results
        halLiftover_merge_bed_file = f"{self.que}.{self.ref}.merge.bed"
        data = pd.read_csv(halLiftover_bed_file, sep="\t", header=None, names=['Chromosome', 'Start', "End", "name"])
        pr_data = pr.PyRanges(data)
        b_pr = pr_data.merge(by="name", slack=self.max_gap)
        df = b_pr.df
        df.to_csv(halLiftover_merge_bed_file, sep="\t", header=False, index=False)
        del data, pr_data, df

        # data = pd.read_csv(halLiftover_bed_file, sep="\t", header=None, names=['chrID', 'start', "end", "que_cns_id"])
        # grouped = data.groupby('que_cns_id', sort=False)
        # data_li = []
        # for group_name, group_data in grouped:
        #     MyBedTooler = MyBedTool(group_data)
        #     outdata = MyBedTooler.merge_bed(self.max_gap)
        #     outdata["que_cns_id"] = group_name
        #     data_li.append(outdata)
        # df = pd.concat(data_li, ignore_index=True, axis=0)
        # df.to_csv(halLiftover_merge_bed_file, sep="\t", header=False, index=False)

        # 3. The area mapped by queCNS to ref overlaps with the CNS of ref
        a_data = pd.read_csv(self.refCNSsBedFile, sep="\t", header=None, names=["Chromosome", "Start", "End", "name"])
        a_pr = pr.PyRanges(a_data)
        # df.columns = ["Chromosome", "Start", "End", "name"]
        # b_pr = pr.PyRanges(df)
        # del df
        intersect_pr = a_pr.join(b_pr, report_overlap=True)
        intersect_df = intersect_pr.df
        intersect_df.insert(4, "Chromosome_Copy", intersect_df["Chromosome"])
        # 4. Generate CNS map file for halLiftover
        halLiftoverAnchorsFile = f"{self.que}.{self.ref}.halLiftover.anchors"
        intersect_df.columns = ["ref_cns_chrID", "ref_cns_start", "ref_cns_end", "ref_cns_ID",
                                "ref_map_chrID", "ref_map_start", "ref_map_end", "que_cns_ID", "ref_map_cns_overlap_length"]
        intersect_df['ratio1'] = intersect_df['ref_map_cns_overlap_length'] / (intersect_df['ref_cns_end'] - intersect_df['ref_cns_start'])
        intersect_df['ratio2'] = intersect_df['ref_map_cns_overlap_length'] / (intersect_df['ref_map_end'] - intersect_df['ref_map_start'])
        intersect_df['max_ratio'] = intersect_df[['ratio1', 'ratio2']].max(axis=1)
        result = intersect_df[intersect_df['max_ratio'] >= self.overlap_rate]
        result = result[["que_cns_ID", "ref_cns_ID"]]
        result = result.drop_duplicates()
        result.to_csv(halLiftoverAnchorsFile, sep="\t", header=False, index=False)

    @staticmethod
    def bed2fa(bedFile, fastaFile, outFastaFile):
        """  Same result as `bedtools getfasta -fi fastaFile -bed bedFile -fo outFastaFile` """
        # Create dictionary to store sequences
        sequences = {}
        # Read the genome file and store the sequence of each chromosome in a dictionary
        with open(fastaFile, "r") as genome_handle:
            for record in SeqIO.parse(genome_handle, "fasta"):
                sequences[record.id] = record.seq
        # Read BED files and extract sequences
        content = ""
        with open(bedFile, "r") as bed_handle:
            for line in bed_handle:
                line_li = line.strip().split("\t")
                chrom, start, end = line_li[0], line_li[1], line_li[2]
                start = int(start)
                end = int(end)
                sequence = sequences[chrom][start:end]
                content += f">{chrom}:{start}-{end}\n{sequence}\n"
        BasePan.write_to(outFastaFile, content)

    def Run_blastn(self):
        queCNSsFasta = f"{self.que}.CNS.fa"
        self.bed2fa(self.queCNSsBedFile, self.queGenomeFile, queCNSsFasta)
        refCNSsFasta = f"{self.ref}.CNS.fa"
        self.bed2fa(self.refCNSsBedFile, self.refGenomeFile, refCNSsFasta)
        cmd = (f"{self.makeblastdb_cmd} -in {refCNSsFasta} -dbtype nucl -out {self.ref} -parse_seqids "
               f"> makeblastdb.log 2>&1")
        self.cmd_linux(cmd)
        # -evalue 1e-6 evalue cannot be set too low, otherwise many maps will be lost.
        cmd = (f"{self.blastn_cmd} -num_threads {self.blastn_threads} "
               f"-query {queCNSsFasta} "
               f"-out {self.que}.{self.ref}.blastn.fmt6.txt "
               f"-db {self.ref} "
               f"-evalue {self.blastn_evalue} "
               f"-task blastn-short "
               f"-outfmt '6 qseqid sseqid pident length qcovs qcovhsp qcovus mismatch gapopen "
               f"qstart qend sstart send evalue bitscore' > blastn.log 2>&1")
        self.cmd_linux(cmd)

    @staticmethod
    def bed2bed(CNSBedFile, outFile):
        data = pd.read_csv(CNSBedFile, sep="\t", header=None, comment='#', skip_blank_lines=True)
        data[4] = "0"
        data[5] = "+"
        data.to_csv(outFile, sep="\t", header=False, index=False)

    def Run_JCVI(self, fm6_blastn_file):
        que_re_CNS_bed, ref_re_CNS_bed = f"{self.que}.bed", f"{self.ref}.bed"
        blastn_data = pd.read_csv(fm6_blastn_file, sep="\t", header=None)
        blastn_data = blastn_data.drop(columns=[4, 5, 6])
        blastn_data.to_csv(f"{self.que}.{self.ref}.last", sep="\t", header=False, index=False)

        # Convert CNSs bed to jcvi input format
        self.bed2bed(self.queCNSsBedFile, que_re_CNS_bed)
        self.bed2bed(self.refCNSsBedFile, ref_re_CNS_bed)

        cmd = "python -m jcvi.compara.catalog ortholog --no_strip_names {} {} " \
              "--cpus={} > jcvi.compara.catalog.ortholog.log 2>&1".format(self.que, self.ref, self.jcvi_threads)
        self.cmd_linux(cmd)

        cmd = (f"python -m jcvi.compara.synteny screen --simple {self.que}.{self.ref}.anchors {self.que}.{self.ref}.new "
               f"> jcvi.compara.synteny.screen.log 2>&1")
        self.cmd_linux(cmd)

    def filter_map(self, halLiftResult, RawMapFile, outResult, checkOverlap=False):
        """
        This function uses halLiftResult to filter the results of blastn and jcvi
        Args:
            halLiftResult: The result of halLiftover after processing
            RawMapFile: blastn or jcvi results
            outResult: Output results file
            checkOverlap: The default is False, RawMapFile only takes the first two columns.
                          If True, RawMapFile must be the result of blastn.
        Returns:
        """
        def getCnsLength(CNSname):
            pattern = r'^(.*?):(\d+)-(\d+)'
            match = re.match(pattern, CNSname)
            if match:
                chrID, start, end = match.groups()
                return int(end) - int(start)
            else:
                raise Exception(f"The format of {CNSname} is incorrect. It must be chrID:start-end. "
                                f"':' cannot appear in chrID")
        halLift_data = pd.read_csv(halLiftResult, sep="\t", header=None, comment='#', skip_blank_lines=True)
        RawMap_data = pd.read_csv(RawMapFile, sep="\t", header=None, comment='#', skip_blank_lines=True)
        data = pd.merge(halLift_data, RawMap_data, on=[0, 1], how='inner')
        if checkOverlap:
            data["len_0"] = data[0].apply(getCnsLength)
            data["len_1"] = data[1].apply(getCnsLength)
            data = data.loc[(data[3] / data[['len_0', 'len_1']].min(axis=1)) >= self.blastn_overlap_rate]
            data = data.drop(columns=['len_0', 'len_1'])
        data = data.drop_duplicates()
        df = data[[0, 1]]
        df.to_csv(outResult, sep="\t", header=False, index=False)
        data.to_csv(f"{outResult}.OriginalFmt.txt", sep="\t", header=False, index=False)

    def halLiftoverFilter_blastn(self, workDir):
        halLiftResult = os.path.join(workDir, "halLiftoverDir", f"{self.que}.{self.ref}",
                                     f"{self.que}.{self.ref}.halLiftover.anchors")
        RawMapFile = os.path.join(workDir, "blastnDir", f"{self.que}.{self.ref}",
                                  f"{self.que}.{self.ref}.blastn.fmt6.txt")
        outResult = os.path.join(workDir, "blastnDir", f"{self.que}.{self.ref}",
                                 f"{self.que}.{self.ref}.blastn.halLiftoverFilter.anchors")
        self.filter_map(halLiftResult, RawMapFile, outResult, checkOverlap=True)

    def halLiftoverFilter_JCVI(self, workDir):
        halLiftResult = os.path.join(workDir, "halLiftoverDir", f"{self.que}.{self.ref}",
                                     f"{self.que}.{self.ref}.halLiftover.anchors")
        RawMapFile = os.path.join(workDir, "JCVIDir", f"{self.que}.{self.ref}",
                                  f"{self.que}.{self.ref}.lifted.anchors")
        outResult = os.path.join(workDir, "JCVIDir", f"{self.que}.{self.ref}",
                                 f"{self.que}.{self.ref}.halLiftoverFilter.lifted.anchors")
        self.filter_map(halLiftResult, RawMapFile, outResult)
        highQualityMapFile = os.path.join(workDir, "JCVIDir", f"{self.que}.{self.ref}",
                                          f"{self.que}.{self.ref}.anchors")
        highQualityOutResult = os.path.join(workDir, "JCVIDir", f"{self.que}.{self.ref}",
                                            f"{self.que}.{self.ref}.halLiftoverFilter.anchors")
        self.filter_map(halLiftResult, highQualityMapFile, highQualityOutResult)

    @staticmethod
    def parse_TowWey_fm6(in_file, out_file):
        data = pd.read_csv(in_file, sep="\t", header=None, comment='#', skip_blank_lines=True)
        columns = list(data.columns)
        columns[0], columns[1] = columns[1], columns[0]
        data = data[columns]
        data.to_csv(out_file, sep="\t", header=False, index=False)

    def get_TowWey_blastn(self):
        filter_raw_file = os.path.join(self.workDir, "blastnDir", f"{self.que}.{self.ref}",
                                       f"{self.que}.{self.ref}.blastn.halLiftoverFilter.anchors.OriginalFmt.txt")
        filter_anchors_file = os.path.join(self.workDir, "blastnDir", f"{self.que}.{self.ref}",
                                           f"{self.que}.{self.ref}.blastn.halLiftoverFilter.anchors")
        outDir = os.path.join(self.workDir, "blastnDir", f"{self.ref}.{self.que}")
        self.pymkdir(outDir)
        out_raw_file = os.path.join(outDir, f"{self.ref}.{self.que}.blastn.halLiftoverFilter.anchors.OriginalFmt.txt")
        out_anchors_file = os.path.join(outDir, f"{self.ref}.{self.que}.blastn.halLiftoverFilter.anchors")
        self.parse_TowWey_fm6(filter_raw_file, out_raw_file)
        self.parse_TowWey_fm6(filter_anchors_file, out_anchors_file)

    def get_TowWey_jcvi(self):
        # filter_raw_file = os.path.join(self.workDir, "JCVIDir", f"{self.que}.{self.ref}",
        #                                f"{self.que}.{self.ref}.halLiftoverFilter.lifted.anchors.OriginalFmt.txt")
        filter_anchors_file = os.path.join(self.workDir, "JCVIDir", f"{self.que}.{self.ref}",
                                           f"{self.que}.{self.ref}.halLiftoverFilter.rescue.lifted.anchors")
        # filter_high_raw_file = os.path.join(self.workDir, "JCVIDir", f"{self.que}.{self.ref}",
        #                                     f"{self.que}.{self.ref}.halLiftoverFilter.anchors.OriginalFmt.txt")
        filter_high_anchors_file = os.path.join(self.workDir, "JCVIDir", f"{self.que}.{self.ref}",
                                                f"{self.que}.{self.ref}.halLiftoverFilter.rescue.anchors")
        outDir = os.path.join(self.workDir, "JCVIDir", f"{self.ref}.{self.que}")
        self.pymkdir(outDir)
        # out_raw_file = os.path.join(outDir, f"{self.ref}.{self.que}.halLiftoverFilter.lifted.anchors.OriginalFmt.txt")
        out_anchors_file = os.path.join(outDir, f"{self.ref}.{self.que}.halLiftoverFilter.rescue.lifted.anchors")
        # out_high_raw_file = os.path.join(outDir, f"{self.ref}.{self.que}.halLiftoverFilter.anchors.OriginalFmt.txt")
        out_high_anchors_file = os.path.join(outDir, f"{self.ref}.{self.que}.halLiftoverFilter.rescue.anchors")
        # self.parse_TowWey_fm6(filter_raw_file, out_raw_file)
        self.parse_TowWey_fm6(filter_anchors_file, out_anchors_file)
        # self.parse_TowWey_fm6(filter_high_raw_file, out_high_raw_file)
        self.parse_TowWey_fm6(filter_high_anchors_file, out_high_anchors_file)

    @staticmethod
    def multi_halLiftover(CNSsAnchor_Instances, outDir):
        # Parallel interface
        os.chdir(outDir)
        CNSsAnchor_Instances.Run_halLiftover()

    @staticmethod
    def multi_RunBlastn(CNSsAnchor_Instances, outDir):
        # Parallel interface
        os.chdir(outDir)
        CNSsAnchor_Instances.Run_blastn()

    @staticmethod
    def multi_halLiftoverFilter_blastn(CNSsAnchor_Instances, outDir, AllWorkDir):
        # Parallel interface
        os.chdir(outDir)
        CNSsAnchor_Instances.halLiftoverFilter_blastn(AllWorkDir)

    @staticmethod
    def multi_JCVI(CNSsAnchor_Instances, outDir, fm6_blastn_file):
        # Parallel interface
        os.chdir(outDir)
        CNSsAnchor_Instances.Run_JCVI(fm6_blastn_file)

    @staticmethod
    def multi_halLiftoverFilter_JCVI(CNSsAnchor_Instances, outDir, AllWorkDir):
        # Parallel interface
        os.chdir(outDir)
        CNSsAnchor_Instances.halLiftoverFilter_JCVI(AllWorkDir)

    @staticmethod
    def multi_get_TowWey(CNSsAnchor_Instances):
        # Parallel interface
        CNSsAnchor_Instances.get_TowWey_blastn()
        CNSsAnchor_Instances.get_TowWey_jcvi()

    def multi_Merge(self, JCVI_filter_File, blastn_filter_File, mergeFile):
        # Parallel interface
        self.merge_files(mergeFile, JCVI_filter_File, blastn_filter_File, uniq=True)
        # cmd = f"cat {JCVI_filter_File} {blastn_filter_File} | sort | uniq > {mergeFile}"
        # self.cmd_linux(cmd)


class MergeSplit:
    """
    First, add the row index to the index_col column of the data frame, then use the columns other than the ref column
    as the index, then group in the same way, and finally retrieve the column according to the initial data frame index.
    This ensures that the final data frame row index obtained for each operation is in the same order.
    """
    def __init__(self, data_a, data_b, merge_col_on):
        self.data1 = copy.deepcopy(data_a)
        self.data2 = data_b
        self.merge_col_on = merge_col_on
        self.check_data()
        self.data1['index_col'] = self.data1.index

    def check_data(self):
        # Determine whether the row index of the data frame is 0,1,2...
        if not (self.data1.index == pd.RangeIndex(start=0, stop=len(self.data1), step=1)).all():
            raise Exception("data1 row index does not match")
        data1_cols = self.data1.columns.tolist()
        data2_cols = self.data2.columns.tolist()
        if self.merge_col_on in data1_cols and self.merge_col_on in data2_cols and len(data2_cols) == 2:
            pass
        else:
            raise Exception("MergeSplit class parameter format error")

    @staticmethod
    def join_with_fillna(series):
        tmp_li = series.fillna('.').tolist()
        if all(element == '.' for element in tmp_li):
            return "."
        else:
            li = [element for element in tmp_li if element != "."]
            return ",".join(li)

    @staticmethod
    def join_merge_col_on(series):
        tmp_li = series.tolist()
        return ",".join(set(tmp_li))

    def get_data1_other_li(self):
        return [i for i in self.data1.columns.tolist() if i != self.merge_col_on]

    def get_data2_other_li(self):
        return [i for i in self.data2.columns.tolist() if i != self.merge_col_on][0]

    def merge_data(self):
        data1_other_cols = self.get_data1_other_li()
        data2_other_col = self.get_data2_other_li()
        data1_split = (self.data1.set_index(data1_other_cols)[self.merge_col_on].str.split(',', expand=True)
                       .stack()
                       .reset_index(level=-1, drop=True)
                       .reset_index(name=self.merge_col_on))
        """
        data1_other_cols = ['B', 'C']
        merge_col_on = 'A'
        
               A      B        C
        0    a,b  aa,bb  aa1,bb1
        1  c,d,e  cc,ff  cc1,ff1
        
        to
        
               B        C  A
        0  aa,bb  aa1,bb1  a
        1  aa,bb  aa1,bb1  b
        2  cc,ff  cc1,ff1  c
        3  cc,ff  cc1,ff1  d
        4  cc,ff  cc1,ff1  e
        """
        merged = pd.merge(data1_split, self.data2, on=self.merge_col_on, how='outer')
        result = merged.groupby(data1_other_cols).agg({
            self.merge_col_on: self.join_merge_col_on,
            data2_other_col: self.join_with_fillna
        }).reset_index()
        result = result.set_index('index_col')
        result = result.reindex(self.data1.index)
        return result[data2_other_col]


class MergeOneCol:
    def __init__(self, data, col_name):
        self.df = copy.deepcopy(data)
        self.df = self.df.reset_index(drop=True)
        self.col_name = col_name

    @staticmethod
    def uniq_key(input_dict):
        while True:
            unique_key = str(uuid.uuid4())
            if unique_key not in input_dict:
                break
        return unique_key

    def cluster_li(self):
        index_li_dict = {}  # Store the index list that needs to be merged
        cell_li_dict = {}  # Store a list of CNS in the index that needs to be merged
        empty_index_li = []
        cell_dict = self.df[self.col_name].to_dict()
        for row_index, cell_txts in cell_dict.items():
            if cell_txts == ".":
                empty_index_li.append(row_index)
                continue
            cell_li = cell_txts.split(",")
            need_merge_index_li = []
            for k, li in cell_li_dict.items():
                if set(cell_li) & set(li):
                    need_merge_index_li.append(k)
            if len(need_merge_index_li) == 0:
                unique_key = self.uniq_key(cell_li_dict)
                cell_li_dict[unique_key] = cell_li
                index_li_dict.setdefault(unique_key, []).append(row_index)
            else:
                merged_cns_li = []  # The merged list.
                merged_cns_li.extend(cell_li)
                merged_index_li = [row_index]
                for i in need_merge_index_li:
                    merged_cns_li.extend(cell_li_dict[i])
                    merged_index_li.extend(index_li_dict[i])
                    del cell_li_dict[i], index_li_dict[i]
                unique_key = self.uniq_key(cell_li_dict)
                cell_li_dict[unique_key] = list(set(merged_cns_li))
                index_li_dict[unique_key] = list(set(merged_index_li))
        return index_li_dict, cell_li_dict, empty_index_li

    @staticmethod
    def merge_col(x):
        if any(val != "." for val in x):
            li = []
            for val in x:
                if val != ".":
                    for j in val.split(","):
                        li.append(j)
            return ','.join(set(li))
        else:
            return "."

    def Merge_on_col(self):
        index_li_dict, _, empty_index_li = self.cluster_li()
        data_li = []
        for i in index_li_dict.keys():
            selected_rows = self.df.loc[index_li_dict[i]]
            # pandas merges each row of the data frame, and each element is connected with ",". If a cell is ".",
            # the cell is ignored, but if the cells of each row in a column are ".", then it is retained. .”
            merged_row = selected_rows.apply(self.merge_col, axis=0)
            new_df = pd.DataFrame([merged_row])
            data_li.append(new_df)
        outDf = pd.concat(data_li, axis=0)
        selected_rows = self.df.loc[empty_index_li]
        outDf = pd.concat([outDf, selected_rows], axis=0)
        return outDf


class RescueCnsBlock(BasePan):
    rate = 0.5

    def __init__(self, cns_block_file, gene_block_file, que_gene_bed_file, ref_gene_bed_file):
        self.cns_block_file = cns_block_file
        self.gene_bed_file = gene_block_file
        self.que_gen_data = self.read_bed_file(que_gene_bed_file)
        self.ref_gen_data = self.read_bed_file(ref_gene_bed_file)

    @staticmethod
    def read_bed_file(bed_file):
        data = pd.read_csv(bed_file, sep='\t', header=None, names=["chrID", "start", "end", "name", "score", "strand"])
        data["ID"] = data.apply(lambda row: "{}:{}-{}".format(row["chrID"], row["start"], row["end"]), axis=1)
        data = data[["name", "ID"]]
        data = data.set_index("name")
        return data

    @staticmethod
    def get_block(row, col_1, col_2):
        chr_1 = re.findall(r"^(.*?):", row[col_1])[0]
        chr_2 = re.findall(r"^(.*?):", row[col_2])[0]
        if chr_1 == chr_2:
            li = [*re.findall(r"[:\-](\d+)", row[col_1]), *re.findall(r"[:\-](\d+)", row[col_2])]
            li = [int(i) for i in li]
            return chr_1, min(li), max(li)
        else:
            raise Exception("chr does not match")

    def parse_gene_block(self):
        gene_block_data = pd.read_csv(self.gene_bed_file, sep="\t", header=None,
                                      names=["que_gene_1", "que_gene_2", "ref_gene_1", "ref_gene_2",
                                             "block_score", "block_strand"])
        gene_block_data = pd.merge(gene_block_data, self.que_gen_data, how="left", left_on="que_gene_1",
                                   right_index=True).rename(columns={'ID': 'que_gene_ID_1'})
        gene_block_data = pd.merge(gene_block_data, self.que_gen_data, how="left", left_on="que_gene_2",
                                   right_index=True).rename(columns={'ID': 'que_gene_ID_2'})

        gene_block_data = pd.merge(gene_block_data, self.ref_gen_data, how="left", left_on="ref_gene_1",
                                   right_index=True).rename(columns={'ID': 'ref_gene_ID_1'})
        gene_block_data = pd.merge(gene_block_data, self.ref_gen_data, how="left", left_on="ref_gene_2",
                                   right_index=True).rename(columns={'ID': 'ref_gene_ID_2'})

        gene_block_data[["que_block_chrID", "que_block_start", "que_block_end"]] = (
            gene_block_data.apply(self.get_block, axis=1, args=("que_gene_ID_1", "que_gene_ID_2")).apply(pd.Series))
        gene_block_data[["ref_block_chrID", "ref_block_start", "ref_block_end"]] = (
            gene_block_data.apply(self.get_block, axis=1, args=("ref_gene_ID_1", "ref_gene_ID_2")).apply(pd.Series))
        gene_block_data = gene_block_data[["que_block_chrID", "que_block_start", "que_block_end",
                                           "ref_block_chrID", "ref_block_start", "ref_block_end"]]
        return gene_block_data

    def parse_cns_block(self):
        cns_block_data = pd.read_csv(self.cns_block_file, sep="\t", header=None,
                                     names=["que_cns_1", "que_cns_2", "ref_cns_1", "ref_cns_2",
                                            "block_score", "block_strand"])
        cns_block_data[["que_block_chrID", "que_block_start", "que_block_end"]] = (
            cns_block_data.apply(self.get_block, axis=1, args=("que_cns_1", "que_cns_2")).apply(pd.Series))
        cns_block_data[["ref_block_chrID", "ref_block_start", "ref_block_end"]] = (
            cns_block_data.apply(self.get_block, axis=1, args=("ref_cns_1", "ref_cns_2")).apply(pd.Series))
        cns_block_data = cns_block_data[["que_block_chrID", "que_block_start", "que_block_end",
                                         "ref_block_chrID", "ref_block_start", "ref_block_end"]]
        return cns_block_data

    def rescue_cns_block(self):
        gene_block_data = self.parse_gene_block()
        cns_block_data = self.parse_cns_block()
        cns_block_data['index'] = range(1, len(cns_block_data) + 1)
        que_cns_block_pr = cns_block_data[["que_block_chrID", "que_block_start", "que_block_end", "index"]]
        que_cns_block_pr.columns = ["Chromosome", "Start", "End", "index"]
        que_cns_block_pr = pr.PyRanges(que_cns_block_pr)

        ref_cns_block_pr = cns_block_data[["ref_block_chrID", "ref_block_start", "ref_block_end", "index"]]
        ref_cns_block_pr.columns = ["Chromosome", "Start", "End", "index"]
        ref_cns_block_pr = pr.PyRanges(ref_cns_block_pr)

        data_li = []
        for index, row in gene_block_data.iterrows():
            row_data = row.to_frame().T

            # que
            que_data = row_data[["que_block_chrID", "que_block_start", "que_block_end"]]
            que_data.columns = ["Chromosome", "Start", "End"]
            que_gene_block_pr = pr.PyRanges(que_data)
            que_overlap_block_pr = que_cns_block_pr.join(que_gene_block_pr, report_overlap=True)
            if que_overlap_block_pr.df.empty:
                index_li1 = []
            else:
                index_li1 = que_overlap_block_pr.df["index"].tolist()

            # ref
            ref_data = row_data[["ref_block_chrID", "ref_block_start", "ref_block_end"]]
            ref_data.columns = ["Chromosome", "Start", "End"]
            ref_gene_block_pr = pr.PyRanges(ref_data)
            ref_overlap_block_pr = ref_cns_block_pr.join(ref_gene_block_pr, report_overlap=True)
            if ref_overlap_block_pr.df.empty:
                index_li2 = []
            else:
                index_li2 = ref_overlap_block_pr.df["index"].tolist()

            intersection = list(set(index_li1) & set(index_li2))
            if len(intersection) > 0:
                que_block_df = que_overlap_block_pr.df[que_overlap_block_pr.df["index"].isin(intersection)]
                que_block_df = que_block_df[["Chromosome", "Start", "End"]]
                que_block_pr = pr.PyRanges(que_block_df)
                que_overlap = que_gene_block_pr.intersect(que_block_pr).sort().merge()
                que_overlap_len = (que_overlap.End - que_overlap.Start).sum()

                ref_block_df = ref_overlap_block_pr.df[ref_overlap_block_pr.df["index"].isin(intersection)]
                ref_block_df = ref_block_df[["Chromosome", "Start", "End"]]
                ref_block_pr = pr.PyRanges(ref_block_df)
                ref_overlap = ref_gene_block_pr.intersect(ref_block_pr).sort().merge()
                ref_overlap_len = (ref_overlap.End - ref_overlap.Start).sum()

                if que_overlap_len / (row['que_block_end'] - row['que_block_start']) > self.rate and ref_overlap_len / (
                        row['ref_block_end'] - row['ref_block_start']) > self.rate:
                    pass
                else:
                    data_li.append(row_data)
            else:
                data_li.append(row_data)
        if len(data_li) > 0:
            result = pd.concat(data_li, axis=0, ignore_index=True)
        else:
            result = pd.DataFrame(columns=["que_block_chrID", "que_block_start", "que_block_end",
                                           "ref_block_chrID", "ref_block_start", "ref_block_end"])
        return result


class RescueCnsAnchor(BasePan):
    flank = 2000

    def __init__(self, block_file, halLiftover_anchors_file):
        self.block_data = self.parse_block_file(block_file)
        self.anchors_data = self.parse_halLiftover_anchors_file(halLiftover_anchors_file)

    def parse_block_file(self, block_file):
        data = pd.read_csv(block_file, sep="\t", header=None,
                           names=["que_block_chrID", "que_block_start", "que_block_end",
                                  "ref_block_chrID", "ref_block_start", "ref_block_end"])
        if len(data) > 0:
            data["que_block_start"] = data["que_block_start"] - self.flank
            data["que_block_end"] = data["que_block_end"] + self.flank
            data["ref_block_start"] = data["ref_block_start"] - self.flank
            data["ref_block_end"] = data["ref_block_end"] + self.flank
        return data

    def parse_halLiftover_anchors_file(self, halLiftover_anchors_file):
        data = pd.read_csv(halLiftover_anchors_file, sep="\t", header=None, names=["que_cns_ID", "ref_cns_ID"])
        data[["que_cns_chrID", "que_cns_start", "que_cns_end"]] = data["que_cns_ID"].apply(
            self.match_CNS_coordinate).apply(pd.Series)
        data[["ref_cns_chrID", "ref_cns_start", "ref_cns_end"]] = data["ref_cns_ID"].apply(
            self.match_CNS_coordinate).apply(pd.Series)
        data = data.drop(columns=["que_cns_ID", "ref_cns_ID"])
        data['index'] = range(1, len(data) + 1)
        return data

    def rescue_cns_anchor(self):
        data_li = []
        if len(self.block_data) > 0:
            for _, row in self.block_data.iterrows():
                row_data = row.to_frame().T

                # que
                que_data = row_data[["que_block_chrID", "que_block_start", "que_block_end"]]
                que_data.columns = ["Chromosome", "Start", "End"]
                que_block_pr = pr.PyRanges(que_data)
                que_cns_pr = self.anchors_data[["que_cns_chrID", "que_cns_start", "que_cns_end", "index"]]
                que_cns_pr.columns = ["Chromosome", "Start", "End", "index"]
                que_cns_pr = pr.PyRanges(que_cns_pr)
                que_overlap_pr = que_cns_pr.join(que_block_pr, report_overlap=True)
                if que_overlap_pr.df.empty:
                    index_li1 = []
                else:
                    index_li1 = que_overlap_pr.df["index"].tolist()

                # ref
                ref_data = row_data[["ref_block_chrID", "ref_block_start", "ref_block_end"]]
                ref_data.columns = ["Chromosome", "Start", "End"]
                ref_block_pr = pr.PyRanges(ref_data)

                ref_cns_pr = self.anchors_data[["ref_cns_chrID", "ref_cns_start", "ref_cns_end", "index"]]
                ref_cns_pr.columns = ["Chromosome", "Start", "End", "index"]
                ref_cns_pr = pr.PyRanges(ref_cns_pr)
                ref_overlap_pr = ref_cns_pr.join(ref_block_pr, report_overlap=True)
                if ref_overlap_pr.df.empty:
                    index_li2 = []
                else:
                    index_li2 = ref_overlap_pr.df["index"].tolist()

                intersection = list(set(index_li1) & set(index_li2))
                if len(intersection) > 0:
                    rescued_data = self.anchors_data[self.anchors_data["index"].isin(intersection)]
                    if len(rescued_data) > 0:
                        data_li.append(rescued_data)
                else:
                    pass
            if len(data_li) > 0:
                result = pd.concat(data_li, axis=0, ignore_index=True)
                result["que"] = result.apply(
                    lambda row: "{}:{}-{}".format(row["que_cns_chrID"], row["que_cns_start"], row["que_cns_end"]), axis=1)
                result["ref"] = result.apply(
                    lambda row: "{}:{}-{}".format(row["ref_cns_chrID"], row["ref_cns_start"], row["ref_cns_end"]), axis=1)
                result = result[["que", "ref"]]
            else:
                result = pd.DataFrame(columns=["que", "ref"])
            return result
        else:
            return pd.DataFrame(columns=["que", "ref"])


def run_rescue_cns_anchor(cnsAnchor, geneWorkDir, configData):
    cns_block_file = os.path.join(cnsAnchor.workDir, "JCVIDir", f"{cnsAnchor.que}.{cnsAnchor.ref}",
                                  f"{cnsAnchor.que}.{cnsAnchor.ref}.simple")
    gene_block_file = os.path.join(geneWorkDir, "JCVIDir", f"{cnsAnchor.que}.{cnsAnchor.ref}",
                                   f"{cnsAnchor.que}.{cnsAnchor.ref}.simple")
    que_gene_bed_file = configData["species"][cnsAnchor.que]["longest_pep_bed"]
    ref_gene_bed_file = configData["species"][cnsAnchor.ref]["longest_pep_bed"]
    RescueCnsBlocker = RescueCnsBlock(cns_block_file, gene_block_file, que_gene_bed_file, ref_gene_bed_file)
    result = RescueCnsBlocker.rescue_cns_block()
    rescue_cns_block_file = os.path.join(cnsAnchor.workDir, "JCVIDir", f"{cnsAnchor.que}.{cnsAnchor.ref}",
                                         f"{cnsAnchor.que}.{cnsAnchor.ref}.rescue_cns_block.bed")
    result.to_csv(rescue_cns_block_file, sep="\t", header=False, index=False)

    halLiftover_anchors_file = os.path.join(cnsAnchor.workDir, "halLiftoverDir", f"{cnsAnchor.que}.{cnsAnchor.ref}",
                                            f"{cnsAnchor.que}.{cnsAnchor.ref}.halLiftover.anchors")
    RescueCnsAnchorer = RescueCnsAnchor(rescue_cns_block_file, halLiftover_anchors_file)
    rescued_data = RescueCnsAnchorer.rescue_cns_anchor()
    rescued_map_file = os.path.join(cnsAnchor.workDir, "JCVIDir", f"{cnsAnchor.que}.{cnsAnchor.ref}",
                                    f"{cnsAnchor.que}.{cnsAnchor.ref}.rescue_cns_map.csv")
    rescued_data.to_csv(rescued_map_file, sep="\t", header=False, index=False)





