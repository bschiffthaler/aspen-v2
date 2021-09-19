#!/bin/bash -l

## error and verbose
set -ex

# modules
module load bioinfo-tools blast

## check for the UPSCb env. var.
if [ -z $UPSCb ]; then
    echo "You need to set the UPSCb environment variable"
    usage
fi

## default args

proj=facility
mail="nicolas.delhomme@slu.se"
in=/mnt/picea/storage/reference/Populus-tremula/v2.2/blast/mRNA/Potra02_transcripts.fasta.gz
out=/mnt/picea/storage/reference/Populus-tremula/v2.2/blast/tmp/mRNA
inx=/mnt/picea/storage/reference/NCBI/20190822/nt
cpu=8

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

## prepare
sbatch -A $proj --mail-user $mail -p all --array=1-400 -c $cpu \
-e $out/Potra02_transcripts.fasta_%a.err -o $out/Potra02_transcripts.fasta_%a.out \
$UPSCb/pipeline/runBlastPlus2.sh -f 5 -p $cpu blastn $in $inx $out 

## then
echo "Combine the results once done: cat Potra02_transcripts.fasta.gz.* > Potra02_transcripts_blast.xml"
