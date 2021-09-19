#!/usr/bin/env Rscript
#' title: "Extract the sequences for B2G"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Set the working dir
setwd("/mnt/picea/storage/reference/Populus-tremula/v2.2")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/storage/reference/Populus-tremula/v2.2")
#' ```

#' Load libraries
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(parallel))

#' Read the mRNA sequence
mRNA <- readDNAStringSet("fasta/Potra02_transcripts.fasta.gz")
prot <- readAAStringSet("fasta/Potra02_proteins.fasta.gz")

#' # Export
#' Prepare the blast input
m_chk <- breakInChunks(length(mRNA),400)
p_chk <- breakInChunks(length(prot),400)

dir.create("blast/mRNA",recursive=TRUE,showWarnings=FALSE)
dev.null <- mclapply(1:length(m_chk),function(i){
  writeXStringSet(mRNA[start(m_chk)[i]:end(m_chk)[i]],
                  file=paste0("blast/mRNA/Potra02_transcripts.fasta.gz.",i))  
},mc.cores=16L)

dir.create("blast/prot",recursive=TRUE,showWarnings=FALSE)
dev.null <- mclapply(1:length(p_chk),function(i){
    writeXStringSet(mRNA[start(p_chk)[i]:end(p_chk)[i]],
                    file=paste0("blast/prot/Potra02_proteins.fasta.gz.",i))  
},mc.cores=16L)


#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```

