#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2024/4/2 15:03
# @Author  : ltan
# @mail    : lei.tan.bio@outlook.com


import re
from Bio import SeqIO
import pandas as pd
import argparse
import logging


class MyParser:
    def __init__(self):
        # self.parser = None
        # self.logger = None
        self.setup_logging()
        self.initialize_parser()
        self.args = self.parser.parse_args()

    def setup_logging(self):
        self.logger = logging.getLogger('MyParser')
        formatter = logging.Formatter('%(asctime)s  %(filename)s[line:%(lineno)d]  %(levelname)s  %(message)s')
        handler = logging.StreamHandler()  # output to the console
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)
        self.logger.setLevel(20)

    def initialize_parser(self):
        self.parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="aaa")
        # self.parser.add_argument('--verbose', help='Common top-level parameter',
        #                     action='store_true', required=False)

        subparsers = self.parser.add_subparsers(help='subcommand', dest='subcommand')

        parent_parser = argparse.ArgumentParser(add_help=False)
        parent_parser.add_argument("-g", "--gff", dest="gffFile", help="The gff files", required=True)
        parent_parser.add_argument("-f", "--fasta", dest="fastaFile", help="The genome fasta files", required=True)
        parent_parser.add_argument("-o", "--prefix", dest="prefix", help="The prefix of output files", required=True)

        parser_gff2fasta = subparsers.add_parser("gff2fasta", parents=[parent_parser], help='gff2fasta')
        parser_gff2fasta.add_argument('-t', "--type", dest="type", choices=['pep', 'cds', 'all'],
                                      help='Choose pep, cds, or all. '
                                           'pep: Extract protein sequence from gff file. '
                                           'cds: Extract cds sequence from gff file. '
                                           'all: Extract protein and cds sequence from gff file.'
                                           'The output files are named {prefix}.pep.fa or {prefix}.cds.fa',
                                      required=True)

        parser_longestGff = subparsers.add_parser("longestGff", parents=[parent_parser], help='longestGff')
        parser_longestGff.add_argument('-t', "--type", dest="type", choices=['cds', 'mRNA'],
                                       help='Choose cds or mRNA. '
                                            'cds: Output transcript annotations for the gene with the longest cds.'
                                            'mRNA: Output transcript annotations for the gene with the longest mRNA.'
                                            'The output files are named {prefix}.longest.cds.gff or {prefix}.longest.mRNA.gff',
                                       required=True)


class GffAnno(MyParser):
    def __init__(self):
        super().__init__()
        self.gffFile = self.args.gffFile
        self.genomeFile = self.args.fastaFile
        self.prefix = self.args.prefix
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
        cds_dict = {}
        for gene_id, mRNA_id_li in self.mRNA_dict.items():
            for mRNA_id in mRNA_id_li:
                cds_dict[mRNA_id] = ""
                data = pd.DataFrame([line.strip().split('\t') for line in self.mRNA_gff_dict[mRNA_id]])
                data.columns = ["chrID", "source", "Type", "start", "end", "score", "strand", "phase", "attributes"]
                cds_data = data[data["Type"] == 'CDS']
                cds_data.loc[:, 'start'] = cds_data['start'].astype(int)
                cds_data.loc[:, 'end'] = cds_data['end'].astype(int)
                # print(cds_data)
                if len(cds_data) == 0:
                    print(f"No CDS in mRNA({mRNA_id})")
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
                print(f"{mRNA_id} not a multiple of three")
            pep_dict[mRNA_id] = str(cds.translate())
        return pep_dict

    @staticmethod
    def split_into_lines(text, width=100):
        lines = []
        for i in range(0, len(text), width):
            lines.append(text[i:i + width])
        return '\n'.join(lines)

    def dict_write_to_fasta(self, file_name, seq_dict):
        txt = ""
        for name, seq in seq_dict.items():
            txt += ">{}\n{}\n".format(name, self.split_into_lines(str(seq)))
        fo = open(file_name, "w")
        fo.write(txt)
        fo.close()

    def output_pep(self):
        output_file = f"{self.prefix}.pep.fa"
        pep_dict = self.extractPep()
        self.dict_write_to_fasta(output_file, pep_dict)

    def output_cds(self):
        output_file = f"{self.prefix}.cds.fa"
        cds_dict = self.extractCds()
        self.dict_write_to_fasta(output_file, cds_dict)

    def run_gff2fasta(self):
        fa_type = self.args.type
        if fa_type == "pep":
            self.output_pep()
        elif fa_type == "cds":
            self.output_cds()
        else:
            self.output_cds()
            self.output_pep()

    def run_longestGff(self):
        longest_type = self.args.type
        txt = ""
        cds_dict = self.extractCds()
        if longest_type == "cds":
            output_gff_file = "{}.longest.cds.gff".format(self.prefix)
            for gene_id, mRNA_id_list in self.mRNA_dict.items():
                tmp_mRNA_dict = {mRNA_id: cds_dict.get(mRNA_id, "") for mRNA_id in mRNA_id_list}
                longest_key = max(tmp_mRNA_dict, key=lambda x: len(tmp_mRNA_dict[x]))
                if len(tmp_mRNA_dict[longest_key]) > 0:
                    txt += self.gene_line_dict[gene_id]
                    txt += self.mRNA_line_dict[longest_key]
                    txt += "".join(self.mRNA_gff_dict[longest_key])
            fo = open(output_gff_file, "w")
            fo.write(txt)
            fo.close()
        elif longest_type == "mRNA":
            output_gff_file = "{}.longest.mRNA.gff".format(self.prefix)
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
                    txt += self.gene_line_dict[gene_id]
                    txt += self.mRNA_line_dict[longest_key]
                    txt += "".join(self.mRNA_gff_dict[longest_key])
            fo = open(output_gff_file, "w")
            fo.write(txt)
            fo.close()

    def run(self):
        if self.args.subcommand == 'gff2fasta':
            self.run_gff2fasta()
        elif self.args.subcommand == 'longestGff':
            self.run_longestGff()
        else:
            print('No subcommand specified')
            raise SystemExit(1)


def main():
    GffAnnoer = GffAnno()
    GffAnnoer.run()


if __name__ == '__main__':
    main()
