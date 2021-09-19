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
in=/mnt/picea/storage/reference/Populus-tremula/v2.2/blast/prot/Potra02_proteins.fasta.gz
out=/mnt/picea/storage/reference/Populus-tremula/v2.2/blast/tmp/prot
inx=/mnt/picea/storage/reference/UniRef90/201908/indices/blast2.6.0+/uniref90
cpu=8

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

## prepare
sbatch -A $proj --mail-user $mail -p all --array=1-400 -c $cpu \
-e $out/Potra02_proteins.fasta_%a.err -o $out/Potra02_proteins.fasta_%a.out \
$UPSCb/pipeline/runBlastPlus2.sh -f 5 -p $cpu blastp $in $inx $out 

## then
echo "Combine the results once done: cat Potra02_proteins.blt.* > Potra02_proteins.xml"
