#'----
#'title: Aspen genome DNA Seq stats
#'author: Bastian Schiffthaler
#'----
library(pander)
library(ShortRead)
library(ggplot2)
library(jsonlite)
library(stringr)
HiSeq.1 <- list(length = 301, type = "Illumina-PE", count = 179989601 * 2,
                tlen = -300, tot = 54176869901 * 2)
HiSeq.2 <- list(length = 301, type = "Illumina-PE", count = 180772624 * 2,
                tlen = -300, tot = 54412559824 * 2)
MiSeq.1 <- list(length = 301, type = "Illumina-PE", count = 28740790 * 2,
                tlen = -300, tot = 8650977790 * 2)
HiSeq.3 <- list(length = 101, type = "Illumina-PE", count = 31308668 * 2,
                tlen = 150, tot = 3162175468 * 2)
HiSeq.4 <- list(length = 101, type = "Illumina-PE", count = 54275444 * 2,
                tlen = 150, tot = 5481819844 * 2)
HiSeq.5 <- list(length = 101, type = "Illumina-PE", count = 101587753 * 2,
                tlen = 300, tot = 10260363053 * 2)
HiSeq.6 <- list(length = 101, type = "Illumina-PE", count = 63533295 * 2,
                tlen = 650, tot = 6416862795 * 2)
roche.454.l <- list(lentgh = NA, type = "454", count = 3317052, tlen = NA,
                    tot = 1727364409)
roche.454.s <- list(length = NA, type = "454", count = 8735661, tlen = NA,
                    tot = 2802021667)
pb.1 <- list(length = NA, type = "PacBio", count = 3581635, tlen = NA,
             tot = 28874072954)
pbf.1 <- list(length = NA, type = "PacBio", count = 131109, tlen = NA,
              tot = 751086642)
pbf.2 <- list(length = NA, type = "PacBio", count = 104220, tlen = NA,
              tot = 826586484)
pbf.3 <- list(length = NA, type = "PacBio", count = 103112, tlen = NA,
              tot = 799030402)
pbf.4 <- list(length = NA, type = "PacBio", count = 104519, tlen = NA,
              tot = 785788626)
pbf.5 <- list(length = NA, type = "PacBio", count = 89210, tlen = NA,
              tot = 623722987)
mp.1 <- list(length = 101, type = "Illumina-MP", count = 105030211 * 2,
             tlen = 3e3, tot = 10608051311 * 2)
mp.2 <- list(length = 101, type = "Illumina-MP", count = 14031890 * 2,
             tlen = 5e3, tot = 1417220890 * 2)
mp.3 <- list(length = 101, type = "Illumina-MP", count = 41611057 * 2,
             tlen = 3e3, tot = 4202716757 * 2)
mp.4 <- list(length = NA, type = "Illumina-MP", count = 70277953 * 2,
             tlen = 5e3, tot = 5753764016)
mp.5 <- list(length = 50, type = "Illumina-MP", count = 27481587 * 2,
             tlen = 1e4, tot = 1374079350 * 2)
fe.1 <- list(length = 151, type = "Lucigen-FE", count = 15602157 * 2,
             tlen = 4e4, tot = 2355925707 * 2)
fe.2 <- list(length = 151, type = "Lucigen-FE", count = 15758362 * 2,
             tlen = 4e4, tot = 2379512662 * 2)
#' ### Table
all.data <- do.call(rbind,lapply(ls(), get))
df <- data.frame(length = unlist(all.data[,1]), type = unlist(all.data[,2]),
                 count = unlist(all.data[,3]), tlen = unlist(all.data[,4]),
                 tot = unlist(all.data[,5]),
                 cov = unlist(all.data[,5]) / 479000000)
pander(df)
pander(sum(df$cov))

#' ### PacBio stats:
pb_files <- dir(paste("/mnt/picea/storage/data/aspseq/nstreet",
                      "/AspenGenomes/DNA-Seq/asp201/PacBio/", sep = ""),
                recursive = TRUE, pattern = "fastq$|fastq.gz$",
                full.names = TRUE)

PBQC <- function(title, dna)
{
  show(ggplot(data.frame(x = width(dna)), aes(x = x)) +
         geom_histogram(binwidth=200) +
         ggtitle(title) + xlab("length"))
  abc <- alphabetFrequency(sread(dna))
  show(ggplot(data.frame(x = (abc[,c("G")] + abc[,"C"]) / rowSums(abc) * 100 ),
              aes(x = x)) + geom_histogram(binwidth = 1) +
         ggtitle(title) + xlab("%GC"))
}
for(f in pb_files)
{
  PBQC(basename(f), readFastq(f))
}

#' ### Detailed 454
pyrodata <- lapply(dir("~/454",pattern = "json", full.names = TRUE), fromJSON)
pyrodata <- unlist(pyrodata)

df <- data.frame(Coverage = pyrodata[grep("Total", names(pyrodata))] / 479000000,
                 MeanRL = pyrodata[grep("Readlength.Mean", names(pyrodata))],
                 MedianRL = pyrodata[grep("Readlength.Median", names(pyrodata))],
                 Count = pyrodata[grep("Readlength.Count", names(pyrodata))],
                 MeanQ = pyrodata[grep("Quality.Mean", names(pyrodata))] - 33,
                 MedianQ = pyrodata[grep("Quality.Median", names(pyrodata))] - 33,
                 MeanGC = pyrodata[grep("GC.Mean", names(pyrodata))],
                 MedianGC = pyrodata[grep("GC.Median", names(pyrodata))])
rownames(df) <- unlist(lapply(strsplit(rownames(df),"\\."), "[[", 1))

ildata <- lapply(dir(paste("/mnt/picea/storage/data/aspseq/nstreet",
                           "/AspenGenomes/DNA-Seq/asp201/summary/",sep = ""),
                     pattern = "json", full.names = TRUE), fromJSON)
ildata <- unlist(ildata)
pb_data <- ildata[grep("pb_", names(ildata))]
ildata <- ildata[grep("pb_", names(ildata), invert = TRUE)]
offset <- ildata[grep("Offset", names(ildata))]

df <- data.frame(Coverage = ildata[grep("Total", names(ildata))] / 479000000,
                 MeanRL = ildata[grep("Readlength.Mean", names(ildata))],
                 MedianRL = ildata[grep("Readlength.Median", names(ildata))],
                 Count = ildata[grep("Readlength.Count", names(ildata))],
                 MeanQ = ildata[grep("Quality.Mean", names(ildata))] - offset,
                 MedianQ = ildata[grep("Quality.Median", names(ildata))] - offset,
                 MeanGC = ildata[grep("GC.Mean", names(ildata))],
                 MedianGC = ildata[grep("GC.Median", names(ildata))])
rownames(df) <- make.unique(rep(c("3KbMP","MiSeq-300","HiSeq-300","10KbMP","5KbMP",
                      "HiSeq-300","3KbMP","5KbMP","PE150","PE150-1",
                      "PE300","PE650","FE1","FE2"), each = 2))
df$Chemistry <- rep(c("old","new","new","old","old","new","old","old","old",
                       "old","old","old","new","new"), each = 2)
ggplot(df, aes(y = MeanGC, fill = Chemistry)) + geom_bar(stat = "identity",
                                                         aes(x = rownames(df))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(df, aes(y = MeanQ, fill = Chemistry)) + geom_bar(stat = "identity",
                                                         aes(x = rownames(df))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

df <- data.frame(Coverage = pb_data[grep("Total", names(pb_data))] / 479000000,
                 MeanRL = pb_data[grep("Readlength.Mean", names(pb_data))],
                 MedianRL = pb_data[grep("Readlength.Median", names(pb_data))],
                 Count = pb_data[grep("Readlength.Count", names(pb_data))],
                 MeanQ = pb_data[grep("Quality.Mean", names(pb_data))],
                 MedianQ = pb_data[grep("Quality.Median", names(pb_data))],
                 MeanGC = pb_data[grep("GC.Mean", names(pb_data))],
                 MedianGC = pb_data[grep("GC.Median", names(pb_data))])
rownames(df) <- unlist(lapply(strsplit(rownames(df),"_fi"), "[[", 1))
pandoc.table(df, style="rmarkdown", split.table = 350)

l <- lapply(dir(paste("/mnt/picea/storage/data/aspseq/nstreet",
                      "/AspenGenomes/DNA-Seq/asp201/summary", sep = ""),
                pattern = "json", full.names = TRUE), fromJSON)
