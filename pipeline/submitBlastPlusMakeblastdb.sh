#!/bin/bash -l

## error and verbose
set -ex

## input
inxDir=/mnt/picea/storage/reference/Populus-tremula/v2.2/indices/BLAST+
fasta=/mnt/picea/storage/reference/Populus-tremula/v2.2/fasta/Potra02_genome.fasta

## create the out dir
if [ ! -d $inxDir ]; then
    mkdir -p $inxDir
fi

## run
sbatch -e $inxDir/Potra02_genome.err -o $inxDir/Potra02_genome.out \
--mail-user="katja.stojkovic@umu.se" $UPSCb/UPSCb-common/pipeline/runBlastPlusMakeblastdb.sh \
$fasta $inxDir
