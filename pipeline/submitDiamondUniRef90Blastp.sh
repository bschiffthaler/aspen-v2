#!/bin/bash

taxonmap=/mnt/picea/storage/reference/Taxonomy/20190825/prot.accession2taxid.gz
taxonnodes=/mnt/picea/storage/reference/Taxonomy/20190825/nodes.dmp
taxonnames=/mnt/picea/storage/reference/Taxonomy/20190825/names.dmp
query=/mnt/picea/storage/reference/Populus-tremula/v2.2/fasta/Potra02_proteins.fasta.gz
inx=/mnt/picea/storage/reference/UniRef90/201908/indices/diamond/uniref90.dmnd
out=/mnt/picea/storage/reference/Populus-tremula/v2.2/diamond
mail=nicolas.delhomme@slu.se

sbatch -t 8:00-00-00 -A facility -e $out/diamond_blastp.err -o $out/diamond_blastp.out --mail-user=$mail \
$UPSCb/pipeline/runDiamond.sh -i $taxonmap -m -n $taxonnodes -s -t $taxonnames blastp $query $inx $out
