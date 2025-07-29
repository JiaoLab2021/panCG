import os
import re
import time
import multiprocessing
import glob
import signal
import subprocess
import pandas as pd
import pyranges as pr
import csv
import shutil
from panCG.lib.base import BasePan
from panCG.lib.base import ParallelScheduler
from panCG.lib.base import GffAnno
from panCG.lib.CnsCall import run_SlideWindow


def calculate_sequence_lengths(fastaFile, outFile):
    content = ""
    Dict_len = {}
    with open(fastaFile, "r") as f:
        for line in f:
            if line.startswith(">"):
                name = re.split(r'[ \t]+', line.strip())[0].replace(">", "")
                Dict_len[name] = 0
            else:
                Dict_len[name] += len(line.strip())
    for name, length in Dict_len.items():
        content += "{}\t{}\n".format(name, length)
    fo = open(outFile, "w")
    fo.write(content)
    fo.close()


def rename_gtf(gtfFile, outGffFile):
    new_gff_content = ""
    with open(gtfFile) as f:
        while True:
            line = f.readline()
            if not line:
                break
            line_li = line.strip().split("\t")
            if line_li[2] == "transcript":
                line_li[2] = "mRNA"
                new_line = "\t".join(line_li) + "\n"
            elif line_li[2] == "exon":
                continue
            elif line_li[2] == "CDS":
                new_line = "\t".join(line_li) + "\n"
            else:
                continue
            new_gff_content += new_line
    f = open(outGffFile, "w")
    f.write(new_gff_content)
    f.close()


def parse_gff(longest_gff_file, out_gff_file=None):
    li = ["mRNA", "exon", "CDS"]

    def get_transcript_id(row):
        if row[2] == "mRNA":
            return re.findall(r"ID=(.*?)(?:;|$)", row[8])[0]
        elif row[2] == "exon":
            return re.findall(r"Parent=(.*?)(?:;|$)", row[8])[0]
        elif row[2] == "CDS":
            return re.findall(r"Parent=(.*?)(?:;|$)", row[8])[0]
    data = pd.read_csv(longest_gff_file, sep="\t", header=None)
    mRNA_data = data[data[2] == "mRNA"].copy()
    mRNA_data["ID"] = mRNA_data[8].apply(lambda x: re.findall(r"ID=(.*?)(?:;|$)", x)[0])
    mRNA_data["Parent"] = mRNA_data[8].apply(lambda x: re.findall(r"Parent=(.*?)(?:;|$)", x)[0])
    mRNA2gene_dict = dict(zip(mRNA_data["ID"], mRNA_data["Parent"]))
    gtf_data = data[data[2].isin(li)].copy()
    gtf_data["mRNA_id"] = gtf_data.apply(get_transcript_id, axis=1)
    gtf_data["gene_id"] = gtf_data["mRNA_id"].apply(lambda x: mRNA2gene_dict.get(x, "-"))
    gtf_data = gtf_data[gtf_data["gene_id"] != "-"]
    gtf_data[8] = gtf_data.apply(lambda row:
                                 'transcript_id "{}"; gene_id "{}"'.format(row["mRNA_id"], row["gene_id"])
                                 if row[2] == "mRNA" else 'transcript_id "{}"; gene_id "{}";'.format(
                                     row["mRNA_id"], row["gene_id"]), axis=1)
    gtf_data = gtf_data.iloc[:, :9]
    gtf_data = gtf_data[gtf_data[2] != "exon"]
    if out_gff_file is not None:
        gtf_data.to_csv(out_gff_file, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE, escapechar="\\")
    return gtf_data


def merge_files(input_files, outputFile):
    with open(outputFile, 'w') as output:
        for file_name in input_files:
            with open(file_name, 'r') as input_file:
                output.write(input_file.read())


def run_command_with_timeout(logger, command, timeout_seconds):
    """ Run an external command and kill it if it is not completed after the time expires """
    process = subprocess.Popen(command, shell=True, preexec_fn=os.setsid)
    try:
        start_time = time.time()
        while process.poll() is None:  # If process.poll() is None, it means that the child process has not ended yet
            time_elapsed = time.time() - start_time
            if time_elapsed > timeout_seconds:
                os.killpg(os.getpgid(process.pid), signal.SIGTERM)
                logger.warning("Process killed due to timeout, cmd: {}".format(command))
                break
            time.sleep(int(timeout_seconds/100))
    except KeyboardInterrupt:
        os.killpg(os.getpgid(process.pid), signal.SIGTERM)
        print("Process killed manually")


def split_wig(wigFile, outDir):
    content = ""
    BasePan.pymkdir(outDir)
    fi = open(wigFile, "r")
    while True:
        line = fi.readline()
        if not line:
            BasePan.write_to(outWig, content)
            break
        if "fixedStep" in line:
            if content != "":
                BasePan.write_to(outWig, content)
            content = ""
            chrID = re.findall(r"chrom=(.+?) ", line)[0]
            start = re.findall(r"start=(.+?) ", line)[0]
            # step = re.findall(r"step=(.+?)", line)[0]
            outWig = os.path.join(outDir, "{}_{}.wig".format(chrID, start))
            content += line
        else:
            content += line
    fi.close()


def split_fasta(input_file, out_dir):
    current_chromosome = None
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if current_chromosome:
                    fo.close()
                current_chromosome = re.split(r'[ \t]+', line.strip())[0].replace(">", "")
                fo = open(os.path.join(out_dir, f'{current_chromosome}.fa'), 'w')
                fo.write(f">{current_chromosome}\n")
            else:
                fo.write(line)
        else:
            fo.close()

def run_callce(logger, config, workDir, Reference, args_dict):
    logger.info("---------------------------------- step 1. read Configuration file ----------------------------------")
    logger.info("args: {}".format(args_dict))
    configData = BasePan.read_yaml(config)
    fastaFile = configData["species"][Reference]["GenomeFile"]
    raw_mafFile = configData["species"][Reference]["mafFile"]
    gffFile = configData["species"][Reference]["gffFile"]
    tree = configData["callCNS"]["tree"]
    speciesList = configData["callCNS"]["speciesList"].split(",")

    timeout_seconds = args_dict["timeout_seconds"]
    estimate_threads = args_dict["estimate_threads"]
    phastCons_threads = args_dict["phastCons_parallel"]

    target_coverage = args_dict["target_coverage"]
    expected_length = args_dict["expected_length"]

    high = args_dict["high"]
    low = args_dict["low"]
    window = args_dict["window"]

    minCol = args_dict["minCol"]
    minRow = args_dict["minRow"]

    mafFilter_cmd = BasePan.check_cmd("mafFilter")
    wigToBigWig_cmd = BasePan.check_cmd("wigToBigWig")
    mafSplit_cmd = BasePan.check_cmd("mafSplit")
    msa_view_cmd = BasePan.check_cmd("msa_view")
    phyloFit_cmd = BasePan.check_cmd("phyloFit")
    msa_split_cmd = BasePan.check_cmd("msa_split")
    phastCons_cmd = BasePan.check_cmd("phastCons")
    phyloBoot_cmd = BasePan.check_cmd("phyloBoot")

    chrList = configData["species"][Reference]["chrList"].split(",")  # chrID used to build the model
    BasePan.pymkdir(workDir)
    BasePan.pymkdir(os.path.join(workDir, "00-genome_split"))
    BasePan.pymkdir(os.path.join(workDir, "01-maf"))
    BasePan.pymkdir(os.path.join(workDir, "02-model"))
    BasePan.pymkdir(os.path.join(workDir, "03-phastCons"))

    logger.info("---------------------------------- step 2. split maf and gff File ----------------------------------")
    os.chdir(workDir)
    logger.info("WorkDir: {}".format(workDir))

    mafFile = os.path.join(workDir, f"{Reference}.filter.maf")
    cmd = f"{mafFilter_cmd} -minCol={minCol} -minRow={minRow} -needComp={Reference} {raw_mafFile} > {mafFile}"
    logger.info(f"runing: {cmd}")
    BasePan.cmd_linux(cmd)

    cmd = f"{mafSplit_cmd} _.bed 01-maf/ {mafFile} -byTarget -useFullSequenceName"
    logger.info(f"runing: {cmd}")
    BasePan.cmd_linux(cmd)
    os.remove(mafFile)
    del mafFile

    longestGffFile = "{}.longest.gff".format(Reference)
    GffAnnoer = GffAnno(logger, gffFile, fastaFile)
    GffAnnoer.longestGff(longestGffFile)
    # outGffFile = "{}.CDS.gff".format(Reference)
    # gtf_data = parse_gff(longestGffFile, outGffFile)
    gtf_data = parse_gff(longestGffFile)
    os.remove(longestGffFile)

    gtf_data[0] = gtf_data[0].apply(lambda x: f"{Reference}.{x}")
    for chrID in chrList:
        out = os.path.join(workDir, "02-model", "{}.CDS.gff".format(chrID))
        df = gtf_data[gtf_data[0] == f"{Reference}.{chrID}"]
        df.to_csv(out, sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE, escapechar="\\")

    logger.info("----------------------------- step 3. Building non-conservative models -----------------------------")
    os.chdir(os.path.join(workDir, "02-model"))
    logger.info("WorkDir: {}".format(os.path.join(workDir, "02-model")))

    BasePan.pymkdir(os.path.join(workDir, "02-model", "logs"))
    BasePan.pymkdir(os.path.join(workDir, "02-model", "CHUNKS"))
    pool = multiprocessing.Pool(phastCons_threads)
    for chrID in chrList:
        cmd = "{} {} --in-format MAF --4d --features {}.CDS.gff > {}_codons.ss"\
            .format(msa_view_cmd, os.path.join(workDir, "01-maf", "{}.maf".format(chrID)), chrID, chrID)
        pool.apply_async(BasePan.cmd_linux, (cmd,))
        time.sleep(0.2)
    pool.close()
    pool.join()

    pool = multiprocessing.Pool(phastCons_threads)
    for chrID in chrList:
        cmd = f"{msa_view_cmd} {chrID}_codons.ss --in-format SS --out-format SS --tuple-size 1 > {chrID}_sites.ss"
        pool.apply_async(BasePan.cmd_linux, (cmd,))
        time.sleep(0.2)
    pool.close()
    pool.join()

    speciesList_str = ",".join(speciesList)
    cmd = f"{msa_view_cmd} --unordered-ss --out-format SS --aggregate {speciesList_str} *_sites.ss > ALL_4d.sites.ss"
    logger.info(f"runing: {cmd}")
    BasePan.cmd_linux(cmd)
    cmd = f"{phyloFit_cmd} --tree \"{tree}\" --msa-format SS --out-root nonconserved_4d ALL_4d.sites.ss " \
          f"> logs/phyloFit.log 2>&1"
    logger.info(f"runing: {cmd}")
    BasePan.cmd_linux(cmd)

    logger.info("-------------------------------- step 4. Estimate conservative model --------------------------------")
    mafFileList = []
    for chrID in chrList:
        mafFileList.append(os.path.join(workDir, "01-maf", f"{chrID}.maf"))

    split_fasta(fastaFile, os.path.join(workDir, "00-genome_split"))
    for mafFile in mafFileList:
        chrID = re.findall(r"^(.+?)\.maf", os.path.basename(mafFile))[0]
        refseq_fa = os.path.join(workDir, "00-genome_split", f"{chrID}.fa")
        cmd = (f"{msa_split_cmd} {mafFile} --in-format MAF --refseq {refseq_fa} --windows 1000000,0 "
               f"--out-root CHUNKS/{chrID} --out-format SS --min-informative 1000 --between-blocks 5000 "
               f"> logs/msa_split.log 2>&1")
        logger.info(f"runing: {cmd}")
        BasePan.cmd_linux(cmd)

    ssFileList = glob.glob(os.path.join(workDir, "02-model", "CHUNKS") + "/*.ss")
    pool = multiprocessing.Pool(estimate_threads)
    for ssFile in ssFileList:
        regionID = re.findall(r"^(.+?)\.ss", os.path.basename(ssFile))[0]
        cmd = f"{phastCons_cmd} --estimate-rho CHUNKS/{regionID}.rho --target-coverage {target_coverage} " \
              f"--expected-length {expected_length} --no-post-probs --msa-format SS " \
              f"--log logs/{regionID}.log {ssFile} nonconserved_4d.mod > logs/{regionID}.nohup.log 2>&1"
        pool.apply_async(run_command_with_timeout, (logger, cmd, timeout_seconds,))
        time.sleep(0.1)
    pool.close()
    pool.join()

    os.chdir(os.path.join(workDir, "02-model", "CHUNKS"))
    logger.info("WorkDir: {}".format(os.path.join(workDir, "02-model", "CHUNKS")))

    logger.info("merge conservative and non-conservative models ...")
    os.chdir(os.path.join(workDir, "02-model", "CHUNKS"))
    logger.info("WorkDir: {}".format(os.path.join(workDir, "02-model", "CHUNKS")))
    cmd = "ls *.rho.cons.mod > cons.txt"
    BasePan.cmd_linux(cmd)
    cmd = f"{phyloBoot_cmd} --read-mods '*cons.txt' --output-average all.ave.cons.mod > phyloBoot.cons.log"
    logger.info(f"runing: {cmd}")
    BasePan.cmd_linux(cmd)

    cmd = "ls *.rho.noncons.mod > noncons.txt"
    BasePan.cmd_linux(cmd)
    cmd = f"{phyloBoot_cmd} --read-mods '*noncons.txt' --output-average all.ave.noncons.mod > phyloBoot.noncons.log"
    logger.info(f"runing: {cmd}")
    BasePan.cmd_linux(cmd)

    logger.info("--------------------------- step 5. phastCons predicts conserved elements ---------------------------")
    os.chdir(os.path.join(workDir, "03-phastCons"))
    logger.info("WorkDir: {}".format(os.path.join(workDir, "03-phastCons")))  
    allMafFileList = glob.glob(os.path.join(workDir, "01-maf") + "/*.maf")
    cmd = "ln -s {} ./".format(os.path.join(workDir, "02-model", "CHUNKS", "all.ave.cons.mod"))
    BasePan.cmd_linux(cmd)
    cmd = "ln -s {} ./".format(os.path.join(workDir, "02-model", "CHUNKS", "all.ave.noncons.mod"))
    BasePan.cmd_linux(cmd)
    
    BasePan.pymkdir(os.path.join(workDir, "03-phastCons", "BED"))
    BasePan.pymkdir(os.path.join(workDir, "03-phastCons", "Wig"))
    BasePan.pymkdir(os.path.join(workDir, "03-phastCons", "logs"))
    pool = multiprocessing.Pool(phastCons_threads)
    for mafFile in allMafFileList:
        chrID = re.findall(r"^(.+?)\.maf", os.path.basename(mafFile))[0]
        cmd = f"{phastCons_cmd} --target-coverage {target_coverage} --expected-length {expected_length} " \
              f"--most-conserved BED/{chrID}.bed --score " \
              f"{mafFile}  all.ave.cons.mod,all.ave.noncons.mod " \
              f"> Wig/{chrID}.wig 2> logs/{chrID}.phastCons.log"
        pool.apply_async(BasePan.cmd_linux, (cmd,))
        time.sleep(0.2)
    pool.close()
    pool.join()

    os.chdir(os.path.join(workDir, "03-phastCons", "Wig"))
    logger.info("WorkDir: {}".format(os.path.join(workDir, "03-phastCons", "Wig")))
    wig_file_li = glob.glob("*.wig")
    all_wig_file = "all.wig"
    with open(all_wig_file, "w") as outfile:
        for wig_file_name in wig_file_li:
            with open(wig_file_name, "r") as infile:
                outfile.write(infile.read())
    falengFile = "{}.chrom.sizes".format(Reference)
    logger.info("create size file ...")
    calculate_sequence_lengths(fastaFile, falengFile)
    cmd = f"{wigToBigWig_cmd} {all_wig_file} {falengFile} {Reference}.all.bw"
    logger.info(f"runing: {cmd}")
    BasePan.cmd_linux(cmd)
    os.remove(falengFile)
    shutil.move(
        os.path.join(workDir, "03-phastCons", "Wig", f"{Reference}.all.bw"),
        os.path.normpath(os.path.join(workDir, "03-phastCons", "Wig", "..", f"{Reference}.all.bw"))
    )

    logger.info("----------------------------- step 6. Sliding window CNS identification -----------------------------")
    wigSplitDir = os.path.join(workDir, "03-phastCons", "Wig", "wigSplitDir")
    split_wig(all_wig_file, wigSplitDir)
    os.remove(all_wig_file)
    BasePan.pymkdir(os.path.join(workDir, "03-phastCons", "Wig", "CEsDir"))
    BasePan.pymkdir(os.path.join(workDir, "03-phastCons", "Wig", "CEsDir", "CEsbedDir"))

    os.chdir(os.path.join(workDir, "03-phastCons", "Wig", "CEsDir"))
    logger.info("WorkDir: {}".format(os.path.join(workDir, "03-phastCons", "Wig", "CEsDir")))
    wigFileList = glob.glob(wigSplitDir + "/*.wig")

    threads = args_dict["threads"]
    total_tasks = len(wigFileList)
    logger.info("total_tasks: {}".format(total_tasks))
    logger.info("threads: {}".format(threads))
    completed_tasks = multiprocessing.Manager().Value('i', 0)
    results = []
    pool = multiprocessing.Pool(processes=threads)
    progress_block = 5000
    parallelScheduler = ParallelScheduler(logger, total_tasks, completed_tasks, progress_block)
    out_bed_file_li = []
    for wigFile in wigFileList:
        name = re.findall(r"^(.+?)\.wig", os.path.basename(wigFile))[0]
        bedFile = os.path.join(workDir, "03-phastCons", "Wig", "CEsDir", "CEsbedDir", f"{name}.bed")
        out_bed_file_li.append(bedFile)
        track_progress = parallelScheduler.make_call_back()
        error_callback = parallelScheduler.make_error_callback()
        result = pool.apply_async(run_SlideWindow, args=(high, low, wigFile, window, bedFile,),
                                  callback=track_progress, error_callback=error_callback)
        results.append(result)
    pool.close()
    pool.join()

    CEsbed_files = glob.glob("CEsbedDir" + "/*.bed")
    merge_files(CEsbed_files, "allCEs.bed")
    allCEs_df = pd.read_csv("allCEs.bed", sep="\t", header=None, names=["Chromosome", "Start", "End"])
    allCEs_pr = pr.PyRanges(allCEs_df)
    allCEs_merge_pr = allCEs_pr.sort().merge(slack=args_dict["merge"])
    allCEs_merge_df = allCEs_merge_pr.df
    allCEs_merge_df.to_csv("allCEs.sort.merge.bed", sep="\t", header=False, index=False)

    # rm tmp file
    for i in wigFileList:
        os.remove(i)
    os.rmdir(wigSplitDir)
    for i in CEsbed_files:
        os.remove(i)
    os.rmdir("CEsbedDir")

    os.remove("allCEs.bed")
        
    for i in glob.glob(os.path.join(workDir, "00-genome_split", "*.fa")):
        os.remove(i)
    os.rmdir(os.path.join(workDir, "00-genome_split"))
    
    for i in glob.glob(os.path.join(workDir, "01-maf", "*.maf")):
        os.remove(i)
    os.rmdir(os.path.join(workDir, "01-maf"))
    
    for i in glob.glob(os.path.join(workDir, "02-model", "*.gff")):
        os.remove(i)
    
    for i in glob.glob(os.path.join(workDir, "02-model", "*.ss")):
        os.remove(i)

    for i in glob.glob(os.path.join(workDir, "02-model", "CHUNKS", "*.ss")):
        os.remove(i)
    
    for i in glob.glob(os.path.join(workDir, "02-model", "CHUNKS", "*.rho.noncons.mod")):
        os.remove(i)

    for i in glob.glob(os.path.join(workDir, "02-model", "CHUNKS", "*.rho.cons.mod")):
        os.remove(i)

    for i in glob.glob(os.path.join(workDir, "03-phastCons", "Wig", "*.wig")):
        os.remove(i)

    for i in glob.glob(os.path.join(workDir, "03-phastCons", "BED", "*.bed")):
        os.remove(i)
    os.rmdir(os.path.join(workDir, "03-phastCons", "BED"))
