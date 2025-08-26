# panCG pipeline
[![PyPI version](https://badge.fury.io/py/panCG.svg)](https://badge.fury.io/py/panCG)

<img src="/figures/panCG.png">

## Dependencies

1. `halLiftover` in [cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/v2.9.3/BIN-INSTALL.md)

2. [phast](https://github.com/CshlSiepelLab/phast)

3. [JCVI](https://github.com/tanghaibao/jcvi)

4. [UCSC](https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/): `mafFilter`, `mafSplit`, `wigToBigWig`
``` shell
wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/mafFilter
wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/mafSplit
wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
```

5. [orthofinder](https://github.com/davidemms/OrthoFinder)

6. [blast](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/)

7. [diamond](https://github.com/bbuchfink/diamond)

## install
Make sure the above dependencies are installed and added to PATH.
``` shell
pip install panCG
panCG -h
```

## usage
``` shell
usage: panCG [-h] [--version]  ...

    an integrative pipeline for family-level super-pangenome analysis across coding and noncoding sequences.

optional arguments:
  -h, --help     show this help message and exit
  --version      show program's version number and exit

Commands:

    callCns      Identification of CNS
    pangene      build gene index
    pancns       build CNS index
    GenePavAsso  Associating gene-PAVs with phenotypes between species
    GLSS         Identification of Gene lineage-specific Synteny networks
    CLSS         Identification of CNS lineage-specific Synteny networks
    CnsGeneLink  According to the relative position relationship between CNS and gene and the maximum number of species supported by CNS index and gene index, CNS index and gene index are linked.
    CnsSyntenyNet
                 Used to construct SyntenyNet for filtered pan-CNS
```

## Input file format requirements
1. The chromosome ID of the genome cannot contain special characters such as `":", "-", ","`, etc., and no other characters except numbers, letters and "_".
2. In the gff annotation file, it is best to only have `gene, mRNA, exon, cds, and utr` information. And gene must contain the `ID` field, and others must contain the `Parent` field.
3. The bed file of gene must be a standard 6-column bed file. `<chrID> <start> <end> <geneID> <score/0> <chain>`.

## Output
### cns calling
| Directory               | File suffix        | Describe                            |
| ----------------------- | ------------------ | ----------------------------------- |
| {Workdir}/03-phastCons/ | {species}.all.bw   | PhastCons Conservative Scoring File |
| {Workdir}/03-phastCons/ | {species}.CNSs.bed | CNS file of {species}               |

### panCNS
| Directory                          | File suffix                                | Describe                                                     |
| ---------------------------------- | ------------------------------------------ | ------------------------------------------------------------ |
| halLiftoverDir                     | {que}.{ref}.bed                            | que’s CNS position on the ref map <ref_map_chrID> <ref_map_start> <ref_map_end> <que_cnsID> |
| halLiftoverDir                     | {que}.{ref}.merge.bed                      | The position of que’s CNS in the ref map (merge if the distance is less than the threshold) <ref_map_chrID> <ref_map_start> <ref_map_end> <que_cnsID> |
| halLiftoverDir                     | {que}.{ref}.halLiftover.anchors            | The correspondence between que's CNS in the hal multiple sequence alignment file and ref's CNS <que_cnsID> <ref_cnsID> |
| halLiftoverDir                     | {que}.{ref}.merge.bw.bed                   | `{que}.{ref}.merge.bed` add averageBwScore and effecve_len   |
| blastnDir                          | .blastn.fmt6.txt                           | Original blastn alignment file                               |
| blastnDir                          | .blastn.halLiftoverFilter.anchors          | blastn anchors after halLiftover filtering                   |
| blastnDir                          | .blastn.halLiftoverFilter.anchors.fmt6.txt | fmt6 format of `.blastn.halLiftoverFilter.anchors` file      |
| JCVIDir                            | .lifted.anchors                            | Recruits additional anchors file output by JCVI              |
| JCVIDir                            | .anchors                                   | High quality anchors file output by JCVI                     |
| JCVIDir                            | .halLiftoverFilter.lifted.anchors          | JCVI recruits additional anchors after halLiftover filtering |
| JCVIDir                            | .halLiftoverFilter.anchors                 | JCVI high quality anchors after halLiftover filtering        |
| {Workdir}/Index/                   | cnsCluster.csv                             | Clustering information of CNS of all species                 |
| {Workdir}/Index/{species}/         | {species}.csv                              | Clustering information of all {species} CNS before           |
| {Workdir}/CEsDir/                  | {species}.recall_cds.bed                   | recall_cds coordinates for each species                      |
| {Workdir}/Ref\_{species}_IndexDir/ | .cnsIndexAssign.csv                        | Results of cnsIndexAssign                                    |
| {Workdir}/Ref\_{species}_IndexDir/ | .cnsIndexMerge.csv                         | Results of cnsIndexMerge                                     |
| {Workdir}/Ref\_{species}_IndexDir/ | .recallCEs.csv                             | The results of recallCEs, which records the results of CE obtained by recall |
| {Workdir}/Ref\_{species}_IndexDir/ | .ReCnsIndexMerge.csv                       | The result after merging `.recallCEs.csv`                    |
| {Workdir}/Ref\_{species}_IndexDir/ | .TripleCnsIndexMerge.csv                   | Classify the CE in `.ReCnsIndexMerge.csv` (cns, cds) and then merge the results |
| {Workdir}/Ref\_{species}_IndexDir/ | .recall.csv                                | The result of merging `.TripleCnsIndexMerge.csv` with cell   |
| {Workdir}/Ref\_{species}_IndexDir/ | .sort.csv                                  | The result of sort `.recall.csv`                             |

### pangene
| Directory                     | File suffix  | Describe           |
| ----------------------------- | ------------ | ------------------ |
| {Workdir}/Ref\_{ref}_IndexDir | .panGene.csv | The result pangene |

The Group column is the homology group identified by orthofinder. 
- OGXXXXXXX.1: Indicates the gene index subdivided in the homology group
- OGXXXXXXX.1.Un: The .Un suffix indicates a set of genes that still exist independently in a single species after CPM.
- OGXXXXXXX.1.tree_1: Indicates the gene index subdivided by gene evolution relationship based on the gene index
- OGXXXXXXX.1.tree_Un: The gene set ending with .tree_Un is a gene set that is not classified using evolutionary relationships.
- UnMapOGXXXXXXX.1: UnMap prefix is the gene that orthofinder has no clustering

## quick start
We provide example data for testing, which can be downloaded at [figshare](https://doi.org/10.6084/m9.figshare.29662034.v1).

### cactus
``` shell
nohup /usr/bin/time -v cactus jobstore species.22way.info.txt Citrus.7ways.test_data.hal \
   --realTimeLogging True \
   --workDir /home/xxx/cactus_dir \
   --maxCores 16 --maxMemory 100G --maxDisk 200G > Citrus.7ways.cactus.log 2>&1 &
   
nohup /usr/bin/time -v cactus-hal2maf jobstore Citrus.7ways.test_data.hal C_sinensis.7ways.maf \
    --refGenome C_sinensis \
    --chunkSize 10000000 \
    --noAncestors \
    --dupeMode single \
    --workDir /home/xxx/cactus_dir > C_sinensis.hal2maf.single.log 2>&1 &
```

### call CNS
```shell
for i in C_sinensis C_limon ponkan C_australasica C_glauca F_hindsii A_buxifolia
do
    /usr/bin/time -v panCG callCns \
        -c /home/ltan/Tmp/01-PanCNSGene_test_data/panCG/Example/CNScalling.config.yaml \
        -w /home/ltan/Tmp/01-PanCNSGene_test_data/01-callcns/${i} \
        -r ${i} > ${i}.callCns.log 2>&1
done
```

### pangene
```shell
nohup /usr/bin/time -v panCG pangene \
    -c /home/ltan/Tmp/01-PanCNSGene_test_data/panCG/Example/panCG.config.yaml \
    -w /home/ltan/Tmp/01-PanCNSGene_test_data/02-pangene \
    -r C_sinensis > pangene.log 2>&1 &
```

### panCNS
```shell
nohup /usr/bin/time -v panCG pancns \
    -c /home/ltan/Tmp/01-PanCNSGene_test_data/panCG/Example/panCG.config.yaml \
    -w /home/ltan/Tmp/01-PanCNSGene_test_data/03-pancns \
    -r C_sinensis \
    -W /home/ltan/Tmp/01-PanCNSGene_test_data/02-pangene \
    > pancns.log 2>&1 &
```

## Citation


