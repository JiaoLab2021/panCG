[![PyPI version](https://badge.fury.io/py/panCG.svg)](https://badge.fury.io/py/panCG)
## 🏛️ Background
PanCG enables the reference-free construction of pangene and pan-conserved noncoding sequences (panCNS) indices and integrative analyses between them.
![Figure1.png30](https://rejo-1323544812.cos.ap-nanjing.myqcloud.com/obsidian/202604082020329.png)

## 🔗 Dependencies

1. `halLiftover` in [cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/v2.9.3/BIN-INSTALL.md)
2. [phast](https://github.com/CshlSiepelLab/phast)
3. [JCVI](https://github.com/tanghaibao/jcvi)
4. [UCSC](https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/): `mafFilter`, `mafSplit`, `wigToBigWig`
``` shell
wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/mafFilter
wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/mafSplit
wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigTo
```
5. [orthofinder](https://github.com/davidemms/OrthoFinder)
6. [blast](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/)
7. [diamond](https://github.com/bbuchfink/diamond)

## 📦 install
确保上述依赖项已安装并添加到 PATH 中。
``` shell
# 直接使用pip安装
pip install panCG
# 测试是否安装成功
panCG -h
```

![image.png30](https://rejo-1323544812.cos.ap-nanjing.myqcloud.com/obsidian/202604081737509.png)

## 🔧 usage
We provide example data for testing, which can be downloaded at [figshare](https://doi.org/10.6084/m9.figshare.29662034.v1).
### 多序列比对
Progressive Cactus [(Armstrong et al. 2020)](https://www.nature.com/articles/s41586-020-2871-y) for multi-genome reference-free alignment
```shell
/usr/bin/time -v cactus jobstore species.22way.info.txt Citrus.7ways.test_data.hal \
   --realTimeLogging True \
   --workDir /home/xxx/cactus_dir \
   --maxCores 16 --maxMemory 100G --maxDisk 200G > Citrus.7ways.cactus.log 2>&1 

/usr/bin/time -v cactus-hal2maf jobstore Citrus.7ways.test_data.hal C_sinensis.7ways.maf \
    --refGenome C_sinensis \
    --chunkSize 10000000 \
    --noAncestors \
    --dupeMode single \
    --workDir /home/xxx/cactus_dir > C_sinensis.hal2maf.single.log 2>&1
```

### callCns
多序列比对之后，我们使用phastCons管道鉴定初始CNS集。
```shell
for i in C_sinensis C_limon ponkan C_australasica C_glauca F_hindsii A_buxifolia
do
    /usr/bin/time -v panCG callcns \
        -c /home/PanCNSGene_test_data/panCG/Example/CNScalling.config.yaml \
        -w /callcns_dir/${i} \
        -r ${i} > ${i}.callCns.log 2>&1
done
```
#### 输入
CNScalling.config.yaml 是输入数据的配置文件，下面是个例子
```yaml
species:
  C_sinensis:
      GenomeFile: /home/01-PanCNSGene_test_data/example_data/C_sinensis/C_sinensis.genome.fasta
      gffFile: /home/01-PanCNSGene_test_data/example_data/C_sinensis/C_sinensis.gff
      mafFile: /home/01-PanCNSGene_test_data/example_data/C_sinensis/C_sinensis.7ways.maf
      chrList: "Chr1"  # If there are multiple chrIDs, use "," split, such as "chr1,chr2,chr3"
...

callCNS:
  speciesList: "C_sinensis,C_limon,ponkan,F_hindsii,C_australasica,C_glauca,A_buxifolia"
  tree: "((((C_sinensis,C_limon),ponkan),((C_australasica,C_glauca),F_hindsii)),A_buxifolia);"
```
其中：
1. 物种名必须和cactus中配置的物种名相同。
2. GenomeFile：基因组文件路径。基因组的染色体ID不能包含“:”、“-”、“,”等特殊字符。换句话说，不能包含除数字、字母和`_`之外的其他字符。
3. gffFile：gff3格式的基因注释文件路径。由于该文件会作为去除与CDS存在overlap的保守序列，所以我们建议如果存在可变剪切应该全部加上。 第三列存在gene、mRNA、exon、CDS、five_prime_UTR（可选）和three_prime_UTR（可选）特征。其中gene的第9列必须包含ID字段，其他的特征必须包含Parent字段。
4. mafFile：多序列比对文件路径。`cactus-hal2maf`输出的maf文件路径。
5. chrList：使用`phyloFit`构建非保守模型的染色体id用逗号分隔。我们建议只保留大于3Mb的染色体id。
6. speciesList：系统发育树上所有物种的列表，用逗号分隔；顺序无关紧要。
7. tree：物种进化关系的拓扑结构。
#### 输出
![image.png30](https://rejo-1323544812.cos.ap-nanjing.myqcloud.com/obsidian/202604081757873.png)

我们关注`/callcns_dir/${species}/03-phastCons/`目录下的结果：
`${species}.all.bw`: 单碱基分辨率下的保守得分(phastCons score)
`allCEs.sort.merge.bed`: 所有的保守序列在基因组的位置
`${species}.CNSs.bed`: 该物种CNS在基因组的位置

### pangene
构建 pangene 并以表格形式显示。列代表每个物种，行代表基因索引。
```shell
/usr/bin/time -v panCG pangene \
    -c /home/PanCNSGene_test_data/panCG/Example/panCG.config.yaml \
    -w /home/PanCNSGene_test_data/02-pangene \
    -r C_sinensis > pangene.log 2>&1
```
注意：这里的`-r`参数并非严格意义上的参考基因组，而是构建pangene的锚点物种。一旦确定之后，接下来的pangene将从该物种的进化距离最近的物种开始，一个接着一个构建。最终的pangene依旧是reference-free的。我们建议选取较晚分化时间的物种作为锚点物种。
#### 输入
`panCG.config.yaml`是输入数据的配置文件，下面是个例子
```yaml
species_tree:
  "(A_buxifolia:0.02787045,((C_australasica:0.0361416,C_glauca:0.0338818)0.388938:0.0152986,((ponkan:0.0125193,F_hindsii:0.0258392)0.187884:0.0133223,(C_sinensis:0.0172991,C_limon:0.0399986)0.19403:0.0104422)0.117647:0.00838269):0.02787045);"

HalFile:
  "/home/01-PanCNSGene_test_data/example_data/Citrus.7ways.test_data.hal"

species:
  A_buxifolia:
    CNS_bed: /home/01-PanCNSGene_test_data/example_data/A_buxifolia/A_buxifolia.CNSs.bed
    GenomeFile: /home/01-PanCNSGene_test_data/example_data/A_buxifolia/A_buxifolia.genome.fasta
    bwFile: /home/01-PanCNSGene_test_data/example_data/A_buxifolia/A_buxifolia.all.bw
    gffFile: /home/01-PanCNSGene_test_data/example_data/A_buxifolia/A_buxifolia.gff
    longest_pep_bed: /home/01-PanCNSGene_test_data/example_data/A_buxifolia/A_buxifolia.longest.bed
    longest_pep_fasta: /home/01-PanCNSGene_test_data/example_data/A_buxifolia/A_buxifolia.longest.pep.fasta
  ...
```
1. species_tree: 带有枝长的进化树，枝长表示物种之间的进化距离。
2. HalFile：cactus比对得到的最初的hal文件
3. CNS_bed：`panCG callcns`得到的CNS文件，即`${species}.CNSs.bed`。
4. GenomeFile：基因组文件路径。基因组的染色体ID不能包含“:”、“-”、“,”等特殊字符。换句话说，不能包含除数字、字母和`_`之外的其他字符。
5. bwFile：`panCG callcns`得到的phastCons保守得分文件，即`${species}.all.bw`。
6. gffFile：gff3格式的基因注释文件路径。由于该文件会作为去除与CDS存在overlap的保守序列，所以我们建议如果存在可变剪切应该全部加上。 第三列存在gene、mRNA、exon、CDS、five_prime_UTR（可选）和three_prime_UTR（可选）特征。其中gene的第9列必须包含ID字段，其他的特征必须包含Parent字段。
7. longest_pep_bed：gffFile中代表性转录本的坐标信息。基因的bed文件必须是标准的6列bed文件，即`<chrID> <start> <end> <gene_ID> <.> <chain>`。
8. longest_pep_fasta: 对应的代表性转录本蛋白序列，名字必须和longest_pep_bed的第4列一致。

#### 输出
![image.png30](https://rejo-1323544812.cos.ap-nanjing.myqcloud.com/obsidian/202604081758375.png)

- 我们关注`/home/PanCNSGene_test_data/02-pangene/Ref_{}_IndexDir`目录：
   `Ref.{}.panGene.final.csv`为最后的pangene。
- `/home/PanCNSGene_test_data/02-pangene/JCVIDir/Synteny.gene.overlap.graph.pkl`是共线性网络文件

### panCNS
基于初始CNS集，整合相似性证据、多序列比对证据和同线性证据来构建panCNS。这也以表格形式呈现。列代表物种，行代表 CNS索引。
```shell
/usr/bin/time -v panCG pancns \
    -c /home/PanCNSGene_test_data/panCG/Example/panCG.config.yaml \
    -w /home/PanCNSGene_test_data/03-pancns \
    -r C_sinensis \
    -W /home/PanCNSGene_test_data/02-pangene \
    > pancns.log 2>&1
```
#### 输入
`-W`参数表示上一步pangene的工作目录，我们发现很小部分实际上的CNS的共线区块可以由于CNS距离过长而被忽略。所以在这里引入gene的共线区块作为一个补充的证据。
`panCG.config.yaml`和pangene一致。
#### 输出
![image.png30](https://rejo-1323544812.cos.ap-nanjing.myqcloud.com/obsidian/202604081758272.png)

- 我们关注`/home/PanCNSGene_test_data/03-pancns/Ref_{}_IndexDir`目录：
 1. `Ref.{}.sort.csv`: 最终的panCNS文件

### CnsSyntenyNet
用于构建 CNS共线性网络
```shell
panCG CnsSyntenyNet \
    -p /home/panCG_test/03-pancns/Ref_C_sinensis_IndexDir/Ref.C_sinensis.sort.csv \
    -c /home/panCG_test/example_data/panCG.config.yaml \
    -w /home/panCG_test/04-CnsSyntenyNet \
    -r C_sinensis \
    -C /home/panCG_test/example_data/panCG.config.yaml \
    -W /home/panCG_test/02-pangene
```
#### 输出
![image.png30](https://rejo-1323544812.cos.ap-nanjing.myqcloud.com/obsidian/202604081759184.png)
1. 由于panCNS中存在`recall_cns`,`recall_cds`,`recall_nonCE`所有这里我们需要重新过滤CNS进行共线性分析。
2. `/home/panCG_test/04-CnsSyntenyNet/JCVIDir/Synteny.CNS.graph.pkl`是CNS共线性网络文件。
3. `new_cns.yaml`是过滤之后的CNS文件。

### GLSS
谱系特异性基因共线性网络的鉴定。
```shell
/usr/bin/time -v panCG GLSS \
    --config /home/Example/geneIndex.config.yaml \
    --workDir /home/workDir \
    --reference C_sinensis \
    --lineage_species_file Cultivated.species.list \
    --output Cultivated-specific.tsv \
    --threads 10

```

### CLSS
谱系特异性CNS共线性网络的识别。
```shell
/usr/bin/time -v panCG CLSS \
    --net_file /home/Synteny.CNS.graph.pkl \
    --lineage_species_file Cultivated.species.list \
    --output Cultivated.CLSS.out.txt \
    --threads 10
```

## 🔍 Output
We have listed some important intermediate files and the final result file for explanation, see the [Markdown](https://github.com/JiaoLab2021/panCG/blob/main/Example/output.md) file.

## 📚 Citation
_Please cite the following article if you use panCG in your research_:
- Lei Tan, Shenchao Zhu, Junli Ye, Wen-Biao Jiao, A new super-pangenome pipeline reveals domestication signatures of conserved noncoding sequences in the orange subfamily, _Molecular Biology and Evolution_, Volume 43, Issue 4, April 2026, msag075, [https://doi.org/10.1093/molbev/msag075](https://doi.org/10.1093/molbev/msag075)

