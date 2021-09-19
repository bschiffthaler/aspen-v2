setwd("/mnt/picea/storage/reference/Populus-tremula/v2.2/diamond/")

# egrep 'Iteration_query-def|Hit_id' Potra02_proteins_Uniref90.xml > prot.tmp
prot <- scan("prot.tmp",what="character",sep="\n")

tx <- scan("tx.tmp",what="character",sep="\n")

p.sel <- grep("Iteration_",prot)

# 73765 proteins
length(p.sel)

t.sel <- grep("Iteration_",tx)

# 73765 transcripts
length(t.sel)

p.id <- prot[p.sel][-which(diff(c(p.sel,length(p.sel)))==1)]

# 67181 protein match
length(p.id)

t.id <- tx[t.sel][-which(diff(c(t.sel,length(t.sel)))==1)]

# 70834 transcript match
length(t.id)

# 72177 are matched
length(union(sub(" p.*","",p.id),sub(" t.*","",t.id)))

# 65838 are common
length(intersect(sub(" p.*","",p.id),sub(" t.*","",t.id)))

# 1343 are unique to the protein
length(setdiff(sub(" p.*","",p.id),sub(" t.*","",t.id)))

# 4996 are unique to the transcript
length(setdiff(sub(" t.*","",t.id),sub(" p.*","",p.id)))

k.id <- t.id[match(setdiff(sub(" t.*","",t.id),sub(" p.*","",p.id)),sub(" t.*","",t.id))]
length(k.id)

save.image("20191002.rda")

library(tidyverse)

p <- TRUE 

csize <- 1e6

outfile="../blast/mRNA/blast/Potra02_transcripts_NCBI_NT_no_UniRef.xml"

file.remove(outfile)

f <- function(data,index){
    message(index)
    message(length(data))
    sel <- grep("Iteration_query-def",data)
    vec <- rep(FALSE,length(data))
    
    pos <- which(data[sel] %in% k.id)
    if(p){
        sel <- c(1,sel)
        pos <- c(1,pos+1)
    }
    
    if(length(pos)>=1){
        vec[unlist(lapply(pos,function(l){sel[l]:ifelse(l==length(sel),length(data),sel[l+1]-1)}))] <- TRUE
    }

    if(length(data)!=csize){
        vec[length(data)-c(3:1)] <- TRUE
    }
    
    if (any(vec)){
        dat <- data[vec]
        sel2 <- grep("Iteration_query-def",dat)
        stopifnot(all(dat[sel2] %in% k.id))
        stopifnot(length(sel2)+p == length(pos))
        
        write(dat,append=TRUE,
              file=outfile)
    }
    p <<- vec[length(data)]
}

# crashes because of an embedded nul character.... Found 4987 so far. So 9 are missing
res <- read_lines_chunked(file="../blast/mRNA/blast/Potra02_transcripts_NCBI_NT.xml",
                   chunk_size=csize,
                   callback=ListCallback$new(f),
                   progress=TRUE)

warning("The output file needs editing to remove the last incomplete <Iteration> block")

data <- read_lines(file="../blast/mRNA/blast/Potra02_transcripts_NCBI_NT.xml",skip=204000000)

sel <- grep("Iteration_query-def",data)
vec <- rep(FALSE,length(data))

pos <- which(data[sel] %in% k.id)

vec[unlist(lapply(pos,function(l){sel[l]:ifelse(l==length(sel),length(data),sel[l+1]-1)}))] <- TRUE

vec[length(data)-c(3:1)] <- TRUE

dat <- data[vec]
sel2 <- grep("Iteration_query-def",dat)
stopifnot(all(dat[sel2] %in% k.id))
stopifnot(length(sel2)+p == length(pos))
    
write(dat,append=TRUE,file=outfile)


