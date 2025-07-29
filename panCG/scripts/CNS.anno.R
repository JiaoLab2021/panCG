suppressMessages(library(optparse))

# Define command line parameters
option_list <- list(
  make_option(c("-b", "--bed"), type = "character", dest = "bed", metavar = "character", help = "Input cns bed file"),
  make_option(c("-g", "--gff"), type = "character", dest = "gff", metavar = "character", help = "GFF annotation file of genes in gene index"),
  make_option(c("-o", "--output"), type = "character", dest = "output", metavar = "character", help = "Output file"),

  make_option(c("-p", "--promoter"), type="integer", dest = "promoter", default=2000, help="The length of promoter [default: %default]"),
  make_option(c("-d", "--downstream"), type="integer", dest = "downstream", default=2000, help="The length of downstream [default: %default]")
)

# Create parsing object
opt_parser <- OptionParser(
  usage = "usage: %prog [options]",
  option_list = option_list,
  add_help_option = TRUE,
  description = "This Script is ..."
)

# Parse command line arguments
opts <- parse_args(opt_parser)

# Check if required parameters are provided
if (is.null(opts$bed) || is.null(opts$gff) || is.null(opts$output) || is.null(opts$promoter) || is.null(opts$downstream)) {
  print_help(opt_parser)
  stop("Error: Missing required argument(s)", call. = FALSE)
}

# run
library(ChIPseeker)
library(GenomicFeatures)
library(tidyverse)
library(rtracklayer)

bed_file <- opts$bed
gff_file <- opts$gff
output_file <- opts$output
promoter_len <- opts$promoter
downstream_len <- opts$downstream

options(ChIPseeker.downstreamDistance = downstream_len)

# 导入 GFF 文件
gff <- rtracklayer::import(gff_file, format = "gff3")

# 将 ID 字段提取为转录本名称
mcols(gff)$Name <- mcols(gff)$ID

# 使用修改后的 GRanges 生成 TxDb
txdb <- makeTxDbFromGRanges(gff)

peaks <- readPeakFile(bed_file)
peakAnno <- annotatePeak(peaks,TxDb=txdb,tssRegion=c(-promoter_len, 500))

df <- as.data.frame(peakAnno)
df <- df %>% mutate(annotation = str_replace_all(annotation, " ", "_"))
write.table(df, file=output_file, sep ="\t", row.names = F, col.names = TRUE, quote = F)


# # 各个region的CNS数量统计
# Freq <- peakAnno@annoStat
# Freq$Num <- round((Freq$Frequency/100) * peakAnno@peakNum)
# df <- data.frame(
#   Variable = specie,
#   value = Freq$Num,
#   taxonomy = Freq$Feature
# )
# write.table(df, out_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


