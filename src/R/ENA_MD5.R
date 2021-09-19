library(here)
library(readr)
md5 <- read_table(here("projects/aspen_v2/doc/ENA/MD5.txt"),col_names=c("MD5","Filename"))

md5$Filename <- sub("^\\.\\/","",md5$Filename)

xml <- scan(here("projects/aspen_v2/doc/ENA/UPSC-0182.Run_noMD5.xml"),what="character",sep="\n")

pos <- sapply(md5$Filename,grep,xml)

xml[pos] <- paste0(sub('"/>$',"",xml[pos]),md5$MD5,'"/>')

write(xml,file=here("projects/aspen_v2/doc/ENA/UPSC-0182.Run.xml"))
