library(Biostrings)
library(ggplot2)
library(scales)
library(chromPlot)
library(readr)
library(dplyr)
library(stringr)
library(GenomicRanges)
library(karyoploteR)
v1 <- readDNAStringSet('/mnt/picea/storage/reference/Populus-tremula/v1.1/fasta/Potra01-genome.fa')
v2 <- readDNAStringSet('/mnt/picea/storage/reference/Populus-tremula/v2.2/fasta/Potra02_genome.fasta')

s <- read_tsv('/mnt/picea/projects/aspseq/asp201/v2/synteny/satsuma_summary.chained.out',
              col_names = c('Qname', 'Qstart', 'Qend', 'Tname', 'Tstart', 'Tend', 'Sim', 'Strand'))

sn <- width(v2)
names(sn) <- str_replace(names(v2), ' ', '_')
gt <- GRanges(seqnames = s$Qname, ranges = IRanges(start = s$Qstart + 1, end = s$Qend),
              strand = s$Strand, seqlengths = sn)
gt <- reduce(gt)
sum(width(gt)) / sum(seqlengths(gt))
gtcr <- gt[str_detect(as.character(seqnames(gt)), '^chr')]
sum(width(gtcr)) / sum(seqlengths(gtcr)[1:19])

pdf('~/pv1synteny.pdf', width = 16, height = 9, pointsize = 24)
kt <- tibble(chr = names(v2[1:19]), start = 1, end = width(v2[1:19]))
pl <- karyoploteR::plotKaryotype(toGRanges(as.data.frame(kt)))
kpPlotRegions(pl, gt[str_detect(as.character(seqnames(gt)), '^chr')], col = "#FFAACC")
dev.off()
