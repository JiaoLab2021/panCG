import os
import re
import shutil
import pandas as pd
import pyranges as pr
from Bio import SeqIO
from panCG.lib.base import BasePan


def get_genome_gap(fastaFile, gapBedFile):
    gapRegion = ""
    with open(fastaFile) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            for match in re.finditer('N+', str(record.seq)):
                # print(record.id, match.start(), match.end(), sep='\t')
                gapRegion += "{}\t{}\t{}\n".format(record.id, match.start(), match.end())
    fo = open(gapBedFile, "w")
    fo.write(gapRegion)
    fo.close()


def getCDS_bed(gff_file, out_file):
    data = pd.read_csv(gff_file, sep="\t", comment="#", skip_blank_lines=True, header=None)
    data.iloc[:, 2] = data.iloc[:, 2].str.upper()
    data = data[data[2] == "CDS"]
    data.iloc[:, 2] = data.iloc[:, 2].str.lower()
    data[3] = data[3] - 1
    data["score"] = 0
    data = data[[0, 3, 4, 2, "score", 6]]
    data.to_csv(out_file, sep="\t", header=False, index=False)


def run_CE_RmCDS_to_CNS(logger, config, workDir, Reference, args_dict):
    logger.info("-------------------------------- step 7. Remove CEs overlap with CDS --------------------------------")
    # read Configuration file
    configData = BasePan.read_yaml(config)
    gffFile = configData["species"][Reference]["gffFile"]
    fastaFile = configData["species"][Reference]["GenomeFile"]
    Gap_overlap_rate = args_dict["gap_rate"]

    os.chdir(os.path.join(workDir, "03-phastCons", "Wig", "CEsDir"))
    logger.info("WorkDir: {}".format(os.path.join(workDir, "03-phastCons", "Wig", "CEsDir")))
    CDSFile = "{}.cds.sort.bed".format(Reference)
    getCDS_bed(gffFile, CDSFile)

    # Remove CEs overlap with CDS to CNS (bedtools subtract)
    allCEs_df = pd.read_csv("allCEs.sort.merge.bed", sep="\t", header=None, names=["Chromosome", "Start", "End"])
    allCEs_pr = pr.PyRanges(allCEs_df)
    cds_df = pd.read_csv(CDSFile, sep="\t", header=None, names=["Chromosome", "Start", "End", "Name", "Score", "Strand"])
    cds_pr = pr.PyRanges(cds_df)
    intersect = allCEs_pr.join(cds_pr, report_overlap=True)
    set1 = set(allCEs_df[["Chromosome", "Start", "End"]].apply(tuple, axis=1))
    if len(intersect) == 0:  # no overlap with cds
        set2 = set()
    else:
        intersect_df = intersect.df
        set2 = set(intersect_df[["Chromosome", "Start", "End"]].apply(tuple, axis=1))
    unique_rows = set1 - set2
    cns1_df = pd.DataFrame(list(unique_rows), columns=["Chromosome", "Start", "End"])
    cns1_pr = pr.PyRanges(cns1_df)
    del set1, set2, unique_rows

    # Remove CNS overlap with Gap
    gapBedFile = "gap.bed"
    get_genome_gap(fastaFile, gapBedFile)
    # bedtools subtract -f {Gap_overlap_rate}
    gap_df = pd.read_csv(gapBedFile, sep="\t", header=None, names=["Chromosome", "Start", "End"])
    gap_pr = pr.PyRanges(gap_df)
    set1 = set(cns1_df[["Chromosome", "Start", "End"]].apply(tuple, axis=1))
    intersect2 = cns1_pr.join(gap_pr, report_overlap=True)
    if len(intersect2) == 0:  # no overlap with gap
        set2 = set()
    else:
        intersect2_df = intersect2.df
        intersect2_df = intersect2_df[intersect2_df["Overlap"] / (intersect2_df["End"] - intersect2_df["Start"]) >= Gap_overlap_rate]
        set2 = set(intersect2_df[["Chromosome", "Start", "End"]].apply(tuple, axis=1))
    unique_rows = set1 - set2
    cns_df = pd.DataFrame(list(unique_rows), columns=["Chromosome", "Start", "End"])
    cns_df["name"] = cns_df.apply(lambda row: "{}:{}-{}".format(row["Chromosome"], row["Start"], row["End"]), axis=1)
    cns_df.to_csv(f"{Reference}.CNSs.bed", sep="\t", header=False, index=False)

    os.remove(CDSFile)
    os.remove(gapBedFile)
    shutil.move(
        os.path.join(workDir, "03-phastCons", "Wig", "CEsDir", "allCEs.sort.merge.bed"),
        os.path.normpath(os.path.join(workDir, "03-phastCons", "Wig", "CEsDir", "..", "..", "allCEs.sort.merge.bed"))
    )
    shutil.move(
        os.path.join(workDir, "03-phastCons", "Wig", "CEsDir", f"{Reference}.CNSs.bed"),
        os.path.normpath(os.path.join(workDir, "03-phastCons", "Wig", "CEsDir", "..", "..", f"{Reference}.CNSs.bed"))
    )
    os.rmdir(os.path.join(workDir, "03-phastCons", "Wig", "CEsDir"))
    os.rmdir(os.path.join(workDir, "03-phastCons", "Wig"))

    logger.info("------------------------------------------ END CNS Calling ------------------------------------------")

