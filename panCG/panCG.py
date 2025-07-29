#!/usr/bin/env python

import sys
import os
import argparse
import logging
from datetime import datetime

from panCG.lib import info
from panCG.bin.callCNS import run_callce
from panCG.bin.CE_RmCDS_to_CNS import run_CE_RmCDS_to_CNS
from panCG.bin.cnsMapMerge import run_cnsMapMerge
from panCG.bin.cnsClustering import run_cnsClustering
from panCG.bin.pancnsAssign import run_pancnsAssign
from panCG.bin.pancnsMerge import run_pancnsMerge
from panCG.bin.cnsRecall import run_cnsRecall
from panCG.bin.pancnsSort import run_pancnsSort

from panCG.bin.geneMapMerge import run_geneMapMerge
from panCG.bin.geneClustering import run_geneClustering
from panCG.bin.pangeneAssign import run_pangeneAssign
from panCG.bin.paralogs import parse_paralogs

from panCG.bin.genePAV import GenePavAsso
from panCG.bin.geneLineageSpeciSyn import glss
from panCG.bin.cnsLineageSpeciSyn import clss
from panCG.bin.CnsGeneLink import run_cns_gene_link
from panCG.bin.CnsSyntenyNet import run_CNS_synteny_net


from panCG.lib.base import TimerDecorator

__version__ = info.__version__
__data__ = datetime.now().strftime("%Y/%m/%d")
__author__ = info.__author__
__email__ = info.__email__
__doc__ = info.__doc__


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass


class MyStart:
    def __init__(self):
        self.logger = self.setup_logging()
        self.parser = argparse.ArgumentParser(prog="panCG", description=__doc__, formatter_class=CustomFormatter, add_help=False)
        self.set_help()
        self.subparsers = self.parser.add_subparsers(title="Commands", dest="command", metavar="")
        self.subparsers.required = True
        self.setup_callcns()
        self.setup_pangene()
        self.setup_pancns()
        self.setup_genePavAsso()
        self.set_GLSS()
        self.set_CLSS()
        self.set_CnsGeneLink()
        self.set_CnsSyntenyNet()

    @staticmethod
    def setup_logging():
        logger = logging.getLogger('MyParser')
        logger.error(f"data: {__data__}")
        logger.error(f"version: {__version__}")
        logger.error(f"author: {__author__}")
        logger.error(f"email: {__email__}\n")
        formatter = logging.Formatter('%(asctime)s\t%(levelname)s\t%(message)s', datefmt='%Y-%m-%d %H:%M:%S')
        handler = logging.StreamHandler()
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.setLevel(20)
        return logger

    def set_help(self):
        self.parser.add_argument("-h", "--help", action="help", help="show this help message and exit")
        self.parser.add_argument("--version", action="version", version=__version__,
                                 help="show program's version number and exit")

    def setup_callcns(self):
        callcns_parser = self.subparsers.add_parser('callcns', help='Identification of CNS',
                                                    formatter_class=argparse.RawDescriptionHelpFormatter)
        # input option
        input_option = callcns_parser.add_argument_group(title="required input")
        input_option.add_argument("-c", "--config", dest="config", required=True,
                                  help="The input configuration file")
        input_option.add_argument("-w", "--workDir", dest="workDir", required=True,
                                  help="Output result directory. If it does not exist, it will be created")
        input_option.add_argument("-r", "--reference", dest="reference", required=True,
                                  help="The reference of index")

        # wig option
        wig_option = callcns_parser.add_argument_group(title="wig parameter")
        wig_option.add_argument("--window", dest="window",
                                type=int, default=5,
                                help="The size of the sliding window [default: %(default)s]")
        wig_option.add_argument("--high", dest="high",
                                type=float, default=0.82,
                                help="The phastCons Score of each base in the window must be higher than this"
                                     "value to start the extension [default: %(default)s]")
        wig_option.add_argument("--low", dest="low",
                                type=float, default=0.54,
                                help="Stop extension when encountering a base with phastCons score less than this "
                                     "value [default: %(default)s]")
        wig_option.add_argument("--merge", dest="merge",
                                type=int, default=10,
                                help="The distance between two adjacent CNSs is smaller than this value and they are "
                                     "merged. Same as `-d` of bedtools merge [default: %(default)s]")

        # phastCons option
        phastCons_option = callcns_parser.add_argument_group(title="phastCons parameter")
        phastCons_option.add_argument("--target_coverage", dest="target_coverage",
                                      type=float, default=0.3,
                                      help="target_coverage parameter [default: %(default)s]")
        phastCons_option.add_argument("--expected_length", dest="expected_length",
                                      type=int, default=45,
                                      help="expected_length parameter [default: %(default)s]")
        phastCons_option.add_argument("--timeout_seconds", dest="timeout_seconds",
                                      type=int, default=10800,
                                      help="The longest model fitting time, which will automatically kill the model "
                                           "if it exceeds this time. [default: %(default)s s]")
        phastCons_option.add_argument("--estimate_threads", dest="estimate_threads",
                                      type=int, default=30,
                                      help="threads for fitting models [default: %(default)s]")
        phastCons_option.add_argument("--phastCons_parallel", dest="phastCons_parallel",
                                      type=int, default=10,
                                      help="phastCons conservative scoring thread, related to chromosome number"
                                           " [default: %(default)s]")

        # mafFilter option
        maf_option = callcns_parser.add_argument_group(title="mafFilter parameter")
        maf_option.add_argument("--minCol", dest="minCol",
                                type=int, default=6,
                                help="The smallest col in multiple sequence alignment [default: %(default)s]")
        maf_option.add_argument("--minRow", dest="minRow",
                                type=int, default=3,
                                help="The smallest row in multiple sequence alignment [default: %(default)s]")

        # cns option
        cns_option = callcns_parser.add_argument_group(title="cns parameter")
        cns_option.add_argument("--threads", dest="threads",
                                type=int, default=36,
                                help="number of sliding window threads [default: %(default)s]")
        cns_option.add_argument("--gap_rate", dest="gap_rate",
                                type=float, default=0.3,
                                help="The number of genomic gaps exceeding this ratio will be removed "
                                     "[default: %(default)s]")

    def setup_pangene(self):
        pangene_parser = self.subparsers.add_parser('pangene', help='build gene index',
                                                      formatter_class=argparse.RawDescriptionHelpFormatter)
        # input option
        input_option = pangene_parser.add_argument_group(title="required input")
        input_option.add_argument("-c", "--config", dest="config", required=True,
                                  help="The input configuration file")
        input_option.add_argument("-w", "--workDir", dest="workDir", required=True,
                                  help="Output result directory. If it does not exist, it will be created")
        input_option.add_argument("-r", "--reference", dest="reference", required=True,
                                  help="The reference of index")

        # diamond option
        diamond_option = pangene_parser.add_argument_group(title="diamond parameter")
        diamond_option.add_argument("--diamond_threads", dest="diamond_threads", type=int, default=6,
                                    help="number of diamond threads [default: %(default)s]")
        diamond_option.add_argument("--diamond_parallel", dest="diamond_parallel", type=int, default=10,
                                    help="number of running diamonds simultaneously [default: %(default)s]")
        diamond_option.add_argument("--identity", dest="identity",
                                    help="filter identity threshold [default: %(default)s]", type=float, default=50.00)
        diamond_option.add_argument("--evalue", dest="evalue", type=float, default=1e-10,
                                    help="filter evalue threshold [default: %(default)s]")

        # jcvi option
        jcvi_option = pangene_parser.add_argument_group(title="jcvi parameter")
        jcvi_option.add_argument("--jcvi_threads", dest="jcvi_threads", type=int, default=6,
                                 help="number of jcvi threads [default: %(default)s]")
        jcvi_option.add_argument("--jcvi_parallel", dest="jcvi_parallel", type=int, default=10,
                                 help="number of running jcvi simultaneously [default: %(default)s]")

        # OrthoFinder option
        orthofinder_option = pangene_parser.add_argument_group(title="OrthoFinder parameter")
        orthofinder_option.add_argument("--orthofinder_threads", dest="orthofinder_threads",
                                        type=int, default=10,
                                        help="number of threads [default: %(default)s]")

        # assign option
        assign_option = pangene_parser.add_argument_group(title="assign parameter")
        assign_option.add_argument("--k_clique", dest="k_clique", type=int, default=3,
                                   help="The size of k-clique in Clique Percolation Method (CPM) [default: %(default)s]")
        assign_option.add_argument("--assign_threads", dest="assign_threads", type=int, default=6,
                                   help="number of assign threads [default: %(default)s]")
        assign_option.add_argument("--assign_chunk_size", dest="assign_chunk_size", type=int, default=1000,
                                   help="buffer of threads [default: %(default)s]")
        assign_option.add_argument("--single_abstract_cutoff", dest="single_abstract_cutoff", type=int, default=5,
                                   help="Abstract nodes with less than or equal to `single_abstract_cutoff` species "
                                        "need to allocate edges with other abstract nodes [default: %(default)s]")
        assign_option.add_argument("--max_cliques_num", dest="max_cliques_num", type=int, default=100000,
                                   help="Communities with more than `max_cliques_num number` will not be processed and "
                                        "will be classified according to their evolutionary relationships. [default: %(default)s]")

        # paralogs option
        paralogs_option = pangene_parser.add_argument_group(title="paralogs parameter")
        paralogs_option.add_argument("--min_paralogs_spe_num", dest="min_paralogs_spe_num", type=int, default=8,
                                     help="If there are paralogous genes in species greater than or equal to "
                                          "`min_paralogs_spe_num` in a gene index, it is necessary to fine-tune the "
                                          "distribution based on the gene evolutionary relationship. [default: %(default)s]")

    def setup_pancns(self):
        pancns_parser = self.subparsers.add_parser('pancns', help='build CNS index',
                                                     formatter_class=argparse.RawDescriptionHelpFormatter)
        # input option
        input_option = pancns_parser.add_argument_group(title="required input")
        input_option.add_argument("-c", "--config", dest="config", required=True,
                                  help="The input configuration file")
        input_option.add_argument("-w", "--workDir", dest="workDir", required=True,
                                  help="Output result directory. If it does not exist, it will be created")
        input_option.add_argument("-r", "--reference", dest="reference", required=True,
                                  help="The reference of index")
        # input_option.add_argument("-C", "--geneConfig", dest="geneConfig", required=True,
        #                           help="the panGene configuration file")
        input_option.add_argument("-W", "--geneWorkDir", dest="geneWorkDir", required=True,
                                  help="the panGene work directory")

        # halLiftover option
        halLiftover_option = pancns_parser.add_argument_group(title="halLiftover parameter")
        halLiftover_option.add_argument("--halLiftover_parallel", dest="halLiftover_parallel",
                                        type=int, default=36,
                                        help="number of running halLiftover simultaneously [default: %(default)s]")
        halLiftover_option.add_argument("--aver_bw_score_threshold", dest="aver_bw_score_threshold",
                                        type=float, default=0.685,
                                        help="(high + low) / 2 [default: %(default)s]")
        halLiftover_option.add_argument("--halLiftover_rate", dest="halLiftover_rate",
                                        type=float, default=0.7,
                                        help=" [default: %(default)s]")

        # blastn option
        blastn_option = pancns_parser.add_argument_group(title="halLiftover parameter")
        blastn_option.add_argument("--blastn_threads", dest="blastn_threads",
                                   type=int, default=6,
                                   help="number of blastn threads [default: %(default)s]")
        blastn_option.add_argument("--blastn_parallel", dest="blastn_parallel",
                                   type=int, default=10,
                                   help="number of running blastn simultaneously [default: %(default)s]")
        blastn_option.add_argument("--blastn_evalue", dest="blastn_evalue",
                                   type=float, default=0.01,
                                   help="evalue of blastn [default: %(default)s]")

        # jcvi option
        jcvi_option = pancns_parser.add_argument_group(title="jcvi parameter")
        jcvi_option.add_argument("--jcvi_threads", dest="jcvi_threads", type=int, default=6,
                                 help="number of jcvi threads [default: %(default)s]")
        jcvi_option.add_argument("--jcvi_parallel", dest="jcvi_parallel", type=int, default=10,
                                 help="number of running jcvi simultaneously [default: %(default)s]")

        # map option
        map_option = pancns_parser.add_argument_group(title="map parameter")
        map_option.add_argument("--map_gap", dest="map_gap", type=int, default=10,
                                help="same as the `-d` parameter of bedtools merge [default: %(default)s]")

        # cluster option
        cluster_option = pancns_parser.add_argument_group(title="cluster parameter")
        cluster_option.add_argument("--cluster_threads", dest="cluster_threads", type=int, default=25,
                                    help="number of cluster threads, maximum number of species [default: %(default)s]")

        # assign option
        assign_option = pancns_parser.add_argument_group(title="assign parameter")
        assign_option.add_argument("--assign_threads", dest="assign_threads", type=int, default=10,
                                   help="number of assign threads [default: %(default)s]")
        assign_option.add_argument("--assign_chunk_size", dest="assign_chunk_size", type=int,
                                   default=10000,
                                   help="buffer of threads [default: %(default)s]")
        assign_option.add_argument("--assign_random_seed", dest="assign_random_seed", type=int,
                                   default=1998,
                                   help="random seed of assign [default: %(default)s]")

        # merge option
        merge_option = pancns_parser.add_argument_group(title="merge parameter")
        merge_option.add_argument("--merge_threads", dest="merge_threads", type=int, default=10,
                                  help="number of merge threads [default: %(default)s]")
        merge_option.add_argument("--merge_chunk_size", dest="merge_chunk_size",
                                  type=int, default=10000,
                                  help="buffer of threads [default: %(default)s]")
        merge_option.add_argument("--merge_random_seed", dest="merge_random_seed", type=int,
                                  default=1998,
                                  help="random seed of merge [default: %(default)s]")
        merge_option.add_argument("--min_cns_num", dest="min_cns_num", type=int, default=2,
                                  help="For less than min_cns_num, if there is overlap, they are merged. If there is "
                                       "no overlap, the best blastn alignment is used [default: %(default)s]")

        # recall option
        recall_option = pancns_parser.add_argument_group(title="recall parameter")
        recall_option.add_argument("--recall_threads", dest="recall_threads", type=int, default=36,
                                   help="number of recall threads [default: %(default)s]")
        recall_option.add_argument("--recall_effective_rate", dest="recall_effective_rate",
                                   type=float, default=0.7,
                                   help="the minimum value of effective_len / max(Candidate_ce, exist_cns) during "
                                        "recall [default: %(default)s]")
        recall_option.add_argument("--recall_cns_rate", dest="recall_cns_rate",
                                   type=float, default=0.5,
                                   help="the minimum value of "
                                        "overlap(Candidate_ce, exist_cns) / min(Candidate_ce, exist_cns), "
                                        "if it exceeds, use exist_cns instead [default: %(default)s]")
        recall_option.add_argument("--recall_chunk_size", dest="recall_chunk_size",
                                   type=int, default=10000,
                                   help="buffer of threads [default: %(default)s]")
    
    def setup_genePavAsso(self):
        genePavAsso_parser = self.subparsers.add_parser('GenePavAsso', help='Associating gene-PAVs with phenotypes between species',
                                                        formatter_class=argparse.RawDescriptionHelpFormatter)
        # input option
        input_option = genePavAsso_parser.add_argument_group(title="required input")
        input_option.add_argument("--pan_gene_file", dest="pan_gene_file", required=True,
                                  help="The pan-Gene file obtained from the `pangene`")
        input_option.add_argument("--phenotype_file", dest="phenotype_file", required=True,
                                  help="The phenotype file, includes two columns, species name and continuous phenotypic value")
        input_option.add_argument("--out_file", dest="out_file", required=True,
                                  help="The output file")
        
        # Optional option
        optional_option = genePavAsso_parser.add_argument_group(title="Optional parameter")
        optional_option.add_argument("--min_species", dest="min_species", type=int, default=5,
                                     help="Maximum gene missing value allowed [default: %(default)s]")

    def set_GLSS(self):
        glss_parser = self.subparsers.add_parser('GLSS', formatter_class=argparse.RawDescriptionHelpFormatter,
                                                 help='Identification of Gene lineage-specific Synteny networks'
                                                 )
        # input option
        input_option = glss_parser.add_argument_group(title="required input")
        input_option.add_argument("-c", "--config", dest="config", required=True,
                                  help="The configuration file in `pangene`")
        input_option.add_argument("-w", "--workDir", dest="workDir", required=True,
                                  help="The workDir result directory in `pangene`")
        input_option.add_argument("-r", "--reference", dest="reference", required=True,
                                  help="The reference of the index is only used to output a list of species with "
                                       "increasing evolutionary distance compared to this species. "
                                       "This list only affects the order of the output file column names and "
                                       "has no effect on the results.")
        input_option.add_argument("-l", "--lineage_species_file", dest="lineage_species_file", required=True,
                                  help="Each row is a file for one species")
        input_option.add_argument("-o", "--output", dest="output_file", required=True,
                                  help="The output file")
        input_option.add_argument("-t", "--threads", dest="threads", type=int, default=10,
                                  help="The threads [default: %(default)s]")

        # Optional option
        optional_option = glss_parser.add_argument_group(title="Optional parameter")
        optional_option.add_argument("-k", "--k_clique", dest="k_clique", type=int, default=3,
                                     help="The size of k-clique in Clique Percolation Method (CPM) [default: %(default)s]")
        optional_option.add_argument("-m", "--miss_spe_num", dest="miss_spe_num", type=int, default=0,
                                     help="Maximum number of species allowed for missing collinearity in Lineage species species [default: %(default)s]")
        optional_option.add_argument("--no_orthofinder", dest="no_orthofinder",
                                     action='store_false',  # Add a bool parameter a with a default value of True
                                     help="Whether to perform analysis in the group identified by orthofinder [default: %(default)s]")

    def set_CLSS(self):
        clss_parser = self.subparsers.add_parser('CLSS', formatter_class=argparse.RawDescriptionHelpFormatter,
                                                 help='Identification of CNS lineage-specific Synteny networks'
                                                 )
        # input option
        input_option = clss_parser.add_argument_group(title="required input")
        # input_option.add_argument("-c", "--config", dest="config", required=True,
        #                           help="The configuration file in `pancns`")
        input_option.add_argument("-n", "--net_file", dest="net_file", required=True,
                                  help="The net file")
        # input_option.add_argument("-r", "--reference", dest="reference", required=True,
        #                           help="The reference of the index is only used to output a list of species with "
        #                                "increasing evolutionary distance compared to this species. This list only "
        #                                "affects the order of the output file column names and has no effect on the results.")
        input_option.add_argument("-l", "--lineage_species_file", dest="lineage_species_file", required=True,
                                  help="Each row is a file for one species")
        input_option.add_argument("-o", "--output", dest="output_file", required=True,
                                  help="The output file")
        input_option.add_argument("-t", "--threads", dest="threads", type=int, default=10,
                                  help="The threads [default: %(default)s]")

        # Optional option
        optional_option = clss_parser.add_argument_group(title="Optional parameter")
        optional_option.add_argument("-k", "--k_clique", dest="k_clique", type=int, default=3,
                                     help="The size of k-clique in Clique Percolation Method (CPM) [default: %(default)s]")
        optional_option.add_argument("-m", "--miss_spe_num", dest="miss_spe_num", type=int, default=0,
                                     help="Maximum number of species allowed for missing collinearity in Lineage species species [default: %(default)s]")

    def set_CnsGeneLink(self):
        CnsGeneLink_parser = self.subparsers.add_parser('CnsGeneLink',
                                                        formatter_class=argparse.RawDescriptionHelpFormatter,
                                                        help='According to the relative position relationship between '
                                                             'CNS and gene and the maximum number of species supported '
                                                             'by CNS index and gene index, CNS index and gene index '
                                                             'are linked.'
                                                        )
        # input option
        input_option = CnsGeneLink_parser.add_argument_group(title="required input")
        input_option.add_argument("-c", "--cns_index", dest="cns_index_file", required=True,
                                  help="The output pan-CNS file in `pancns`")
        input_option.add_argument("-g", "--gene_index", dest="gene_index_file", required=True,
                                  help="The output pan-Gene file in `pangene`")
        input_option.add_argument("-l", "--link_yaml", dest="link_yaml", required=True,
                                  help="The contig yaml file. Annotation file recording each gene of each species in the `pangene`")
        input_option.add_argument("-o", "--out_dir", dest="out_dir", required=True,
                                  help="The out dir")

        # Optional option
        optional_option = CnsGeneLink_parser.add_argument_group(title="Optional parameter")
        optional_option.add_argument("-t", "--threads", dest="threads", type=int, default=1,
                                     help="Number of parallel Annotation [default: %(default)s]")

    def set_CnsSyntenyNet(self):
        CnsSyntenyNet_parser = self.subparsers.add_parser('CnsSyntenyNet',
                                                          formatter_class=argparse.RawDescriptionHelpFormatter,
                                                          help='Used to construct SyntenyNet for filtered pan-CNS'
                                                          )
        # input option
        input_option = CnsSyntenyNet_parser.add_argument_group(title="required input")
        input_option.add_argument("-p", "--panCns", dest="panCns", required=True,
                                  help="The panCns file")
        input_option.add_argument("-c", "--cnsConfig", dest="cnsConfig", required=True,
                                  help="The input panCns configuration file")
        input_option.add_argument("-w", "--workDir", dest="workDir", required=True,
                                  help="Output result directory")
        input_option.add_argument("-r", "--reference", dest="reference", required=True,
                                  help="The reference of index")
        input_option.add_argument("-C", "--geneConfig", dest="geneConfig", required=True,
                                  help="the panGene configuration file")
        input_option.add_argument("-W", "--geneWorkDir", dest="geneWorkDir", required=True,
                                  help="the panGene work directory")

        # halLiftover option
        halLiftover_option = CnsSyntenyNet_parser.add_argument_group(title="halLiftover parameter")
        halLiftover_option.add_argument("--halLiftover_parallel", dest="halLiftover_parallel",
                                        type=int, default=36,
                                        help="number of running halLiftover simultaneously [default: %(default)s]")
        halLiftover_option.add_argument("--aver_bw_score_threshold", dest="aver_bw_score_threshold",
                                        type=float, default=0.685,
                                        help="(high + low) / 2 [default: %(default)s]")
        halLiftover_option.add_argument("--halLiftover_rate", dest="halLiftover_rate",
                                        type=float, default=0.7,
                                        help=" [default: %(default)s]")

        # blastn option
        blastn_option = CnsSyntenyNet_parser.add_argument_group(title="halLiftover parameter")
        blastn_option.add_argument("--blastn_threads", dest="blastn_threads",
                                   type=int, default=6,
                                   help="number of blastn threads [default: %(default)s]")
        blastn_option.add_argument("--blastn_parallel", dest="blastn_parallel",
                                   type=int, default=10,
                                   help="number of running blastn simultaneously [default: %(default)s]")
        blastn_option.add_argument("--blastn_evalue", dest="blastn_evalue",
                                   type=float, default=0.01,
                                   help="evalue of blastn [default: %(default)s]")

        # jcvi option
        jcvi_option = CnsSyntenyNet_parser.add_argument_group(title="jcvi parameter")
        jcvi_option.add_argument("--jcvi_threads", dest="jcvi_threads", type=int, default=6,
                                 help="number of jcvi threads [default: %(default)s]")
        jcvi_option.add_argument("--jcvi_parallel", dest="jcvi_parallel", type=int, default=10,
                                 help="number of running jcvi simultaneously [default: %(default)s]")

        # map option
        map_option = CnsSyntenyNet_parser.add_argument_group(title="map parameter")
        map_option.add_argument("--map_gap", dest="map_gap", type=int, default=10,
                                help="same as the `-d` parameter of bedtools merge [default: %(default)s]")

    @staticmethod
    def decorated_func(func1, logger, *args, **kwargs):
        timer = TimerDecorator(logger)
        func2 = timer(func1)
        return func2(*args, **kwargs)

    # --------------------------- cmd --------------------------- #
    def run_callcns(self):
        args = self.parser.parse_args()
        args_dict = vars(args)
        config = args.config
        workDir = args.workDir
        reference = args.reference
        # 创建计时器装饰器（传入当前实例的 logger）
        timer = TimerDecorator(self.logger)
        timer(run_callce)(self.logger, config, workDir, reference, args_dict)
        timer(run_CE_RmCDS_to_CNS)(self.logger, config, workDir, reference, args_dict)

    def run_pangene(self):
        args = self.parser.parse_args()
        args_dict = vars(args)
        config = args.config
        workDir = args.workDir
        reference = args.reference
        timer = TimerDecorator(self.logger)
        timer(run_geneMapMerge)(self.logger, config, workDir, reference, args_dict)
        timer(run_geneClustering)(self.logger, config, workDir, reference, args_dict)
        timer(run_pangeneAssign)(self.logger, config, workDir, reference, args_dict)
        timer(parse_paralogs)(self.logger, config, workDir, reference, args_dict)

    def run_pancns(self):
        args = self.parser.parse_args()
        args_dict = vars(args)
        config = args.config
        workDir = args.workDir
        reference = args.reference
        geneWorkDir = args.geneWorkDir
        # geneConfig = args.geneConfig
        timer = TimerDecorator(self.logger)
        timer(run_cnsMapMerge)(self.logger, config, workDir, geneWorkDir, reference, args_dict)
        timer(run_cnsClustering)(self.logger, config, workDir, reference)
        timer(run_pancnsAssign)(self.logger, config, workDir, reference, args_dict)
        timer(run_pancnsMerge)(self.logger, config, workDir, reference, args_dict)
        timer(run_cnsRecall)(self.logger, config, workDir, reference, args_dict)
        timer(run_pancnsSort)(self.logger, config, workDir, reference)
    
    def run_gene_PAV_asso(self):
        args = self.parser.parse_args()
        pan_gene_file = args.pan_gene_file
        phenotype_file = args.phenotype_file
        out_file = args.out_file
        min_species = args.min_species
        GenePavAsso.set_min_species(min_species)
        GenePavAssoer = GenePavAsso(pan_gene_file, phenotype_file, out_file)
        GenePavAssoer.gene_pav_asso()

    def run_GLSS(self):
        args = self.parser.parse_args()
        args_dict = vars(args)
        config = args.config
        workDir = args.workDir
        reference = args.reference
        lineage_species_file = args.lineage_species_file
        output_file = args.output_file
        threads = args.threads
        glss(self.logger, config, workDir, reference, lineage_species_file, output_file, args_dict, threads)

    def run_CLSS(self):
        args = self.parser.parse_args()
        args_dict = vars(args)
        net_file = args.net_file
        lineage_species_file = args.lineage_species_file
        output_file = args.output_file
        threads = args.threads
        clss(self.logger, net_file, lineage_species_file, output_file, args_dict, threads)

    def run_CnsGeneLink(self):
        args = self.parser.parse_args()
        cns_index_file = args.cns_index_file
        gene_index_file = args.gene_index_file
        link_yaml = args.link_yaml
        out_dir = args.out_dir
        threads = args.threads
        run_cns_gene_link(self.logger, cns_index_file, gene_index_file, link_yaml, threads, out_dir)

    def run_CnsSyntenyNet(self):
        args = self.parser.parse_args()
        args_dict = vars(args)
        panCns = args.panCns
        cnsConfig = args.cnsConfig
        workDir = args.workDir
        reference = args.reference
        geneWorkDir = args.geneWorkDir
        geneConfig = args.geneConfig
        run_CNS_synteny_net(self.logger, panCns, cnsConfig, workDir, reference, geneWorkDir, args_dict)

    # --------------------------- cmd --------------------------- #
    def run(self):
        args = self.parser.parse_args()
        if args.command == 'callcns':
            self.run_callcns()
        elif args.command == 'pancns':
            self.run_pancns()
        elif args.command == 'pangene':
            self.run_pangene()
        elif args.command == 'GenePavAsso':
            self.run_gene_PAV_asso()
        elif args.command == 'GLSS':
            self.run_GLSS()
        elif args.command == 'CLSS':
            self.run_CLSS()
        elif args.command == 'CnsGeneLink':
            self.run_CnsGeneLink()
        elif args.command == 'CnsSyntenyNet':
            self.run_CnsSyntenyNet()
        else:
            self.logger.error('No command specified')
            raise SystemExit(1)


def main():
    MyStartClass = MyStart()
    MyStartClass.run()


if __name__ == '__main__':
    main()

