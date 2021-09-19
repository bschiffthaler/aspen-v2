library(Biostrings)
library(ggplot2)
library(ggrepel)
library(scales)
potra_v2_genome <- readDNAStringSet("/mnt/picea/storage/reference/Populus-tremula/v2.2/fasta/Potra02_genome.fasta")
w52_genome <-  readDNAStringSet("/mnt/picea/storage/reference/Populus-tremula/v2.2/fasta/W52-T195_curated_reference-scaffolded.fasta.gz")


genome <- w52_genome
aln <- w52_chr
seqlevels(aln) <- str_replace(seqlevels(aln), "_RaGOO", "")

o1 <- str_extract(seqnames(aln), "^[csS]")
o2 <- as.numeric(str_extract(seqnames(aln), "\\d+"))
o3 <- start(aln)
o4 <- end(aln)

aln <- aln[order(o1, o2, o3, o4)]

# Uncomment to zoom in on a chr
#aln <- aln[seqnames(aln) %in% "chr19"]

names(genome) <- str_replace(names(genome), "_RaGOO", "")
o1 <- str_extract(names(genome), "^[csS]")
o2 <- as.numeric(str_extract(names(genome), "\\d+"))

genome <- genome[order(o1, o2)]

genome <- genome[names(genome) %in% seqnames(aln)]

cs <- c(0, cumsum(width(genome)))
names(cs) <- names(genome)
cs <- cs[names(cs) %in% as.character(seqnames(aln))]

starts <- cs[as.character(seqnames(aln))] + start(aln)

cols <- rep(c("C1", "C2"), length(cs) / 2)
if (length(cols) < length(cs)) {
  cols <- c(cols, "C1")
}
names(cols) <- names(cs)

chromosomes <- tibble(Start = cs, End = cs + width(genome[names(cs)]),
                      Name = names(cs),
                      Col = cols[names(cs)])


vals <- tibble(X = starts, Y = mcols(aln)[["importance"]],
               Name = as.character(seqnames(aln)), Col = cols[as.character(seqnames(aln))])

ymax <- min(vals$Y)
ymax <- ymax - abs(ymax)
ymin <- ymax - (max(vals$Y) - min(vals$Y)) * 0.075
theme_set(theme_bw(base_size = 18))
ggplot() +
  geom_rect(data = chromosomes, aes(xmin = Start, xmax = End, ymin = ymin, ymax = ymax, fill = Col), col = "white") +
  geom_label(data = chromosomes, aes(label = Name, x = (Start + End) / 2, y = (ymin + ymax) / 2), size = 3) +
  geom_point(data = vals, aes(x = X, y = Y, col = Col, alpha = Y / max(Y))) +
  theme_bw(base_size = 18) +
  theme(legend.position = "none") +
  scale_color_grey() +
  scale_fill_grey() +
  ylab("Importance") +
  xlab("") +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  theme(legend.position = "none") 
# Uncomment t zoom in on a range
# + lim(c(9e6, 975e4))
ggsave("~/sex_male.pdf", width = 16, height = 9, dpi = 300)
