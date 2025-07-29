#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2024/6/1 17:10
# @Author  : ltan
# @mail    : lei.tan.bio@outlook.com

import glob
import os
import yaml


def main():
    li = glob.glob("/home/ltan/datas/03-Brassicaceae/01-genomes/*/*.longest.cds.gff")
    yaml_dict = {"species": {}}
    for i in li:
        nameID = os.path.basename(i).replace(".longest.cds.gff", "")
        yaml_dict["species"][nameID] = {}
        longest_pep_fasta_file = os.path.join(os.path.dirname(i), f"{nameID}.longest.pep.fa")
        longest_pep_bed_file = os.path.join(os.path.dirname(i), f"{nameID}.longest.pep.sort.bed")
        if not os.path.exists(longest_pep_fasta_file) or not os.path.exists(longest_pep_bed_file):
            print(longest_pep_fasta_file)
            raise Exception(f"No longest peptide fasta")

        yaml_dict["species"][nameID]["longest_pep_fasta"] = longest_pep_fasta_file
        yaml_dict["species"][nameID]["longest_pep_bed"] = longest_pep_bed_file

    with open('/home/ltan/datas/03-Brassicaceae/12-geneIndex_20240601/panCG-main_20240601/geneIndex.cnsIndex.config.yaml',
              'w') as file:
        yaml.dump(yaml_dict, file, default_flow_style=False, allow_unicode=True)

    print("The dictionary was successfully written to the YAML file")


if __name__ == '__main__':
    main()
