library(Biostrings)
  pb <- readDNAStringSet("/mnt/picea/storage/data/aspseq/nstreet/AspenGenomes/DNA-Seq/asp201/PacBio/genomic/pb_158/filtered_subreads/pb_158_filtered_subreads.fasta")
f <- lapply(strsplit(names(pb), "/"), "[[", 1)
pbs <- split(pb,unlist(f ))

for(f in names(pbs)){
  writeXStringSet(pbs[[f]],paste("/mnt/picea/projects/aspseq/asp201/processed-data/pb_raw/",
                                 f,".fasta",sep=""))
}