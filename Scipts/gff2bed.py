#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2024/6/1 16:06
# @Author  : ltan
# @mail    : lei.tan.bio@outlook.com


import re
import logging
import argparse
import pandas as pd


class MyParser:
    def __init__(self):
        self.logger = self.setup_logging()
        self.parser = self.initialize_parser()
        self.args = self.parser.parse_args()

    @staticmethod
    def setup_logging():
        logger = logging.getLogger('MyParser')
        formatter = logging.Formatter('%(asctime)s  %(filename)s[line:%(lineno)d]  %(levelname)s  %(message)s')
        handler = logging.StreamHandler()  # output to the console
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.setLevel(20)
        return logger

    @staticmethod
    def initialize_parser():
        parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="aaa")

        parser.add_argument("-f", "--gff", dest="gff",
                            help="The gff file", required=True)
        parser.add_argument("-t", "--type", dest="type",
                            help="The value of the third column of gff", required=True)
        parser.add_argument("-k", "--key", dest="key",
                            help="The key value of the ninth column of gff", required=True)
        parser.add_argument("-o", "--output", dest="output",
                            help="The output", required=True)
        return parser


class Gff2bed(MyParser):
    def __init__(self):
        super().__init__()
        self.gff = self.args.gff
        self.type = self.args.type
        self.key = self.args.key
        self.output = self.args.output

    def run_gff2bed(self):
        data = pd.read_csv(self.gff, sep="\t", header=None, comment='#', skip_blank_lines=True)
        data.columns = ["chrID", "source", "Type", "start", "end", "score", "strand", "phase", "attributes"]
        data = data[data["Type"] == self.type]
        data["nameID"] = data.apply(lambda row: re.findall(r"{}=([^;$]+)".format(self.key), row["attributes"])[0], axis=1)
        data['score'] = 0
        data['start'] = data['start'] - 1
        results = data[["chrID", "start", "end", "nameID", "score", "strand"]]
        results.to_csv(self.output, sep="\t", header=False, index=False)


def main():
    Gff2beder = Gff2bed()
    Gff2beder.run_gff2bed()


if __name__ == '__main__':
    main()
