/home/ltan/biosoft/Cactus/cactus-bin-v2.5.1/bin/halLiftover /home/ltan/Tmp/PanCNSgene_testFile/Test.chr1.6-way.hal Citrus_sinensis <(echo -e "chr1\t9994867\t9994954") Citrus_grandis test1.bed
/home/ltan/.local/bin/bwtool extract bed <(bedtools sort -i test1.bed) /home/ltan/Tmp/PanCNSgene_test_v1.0.1/01-callCNS/Citrus_grandis/03-phastCons/Wig/Citrus_grandis.all.bw test2.bed

