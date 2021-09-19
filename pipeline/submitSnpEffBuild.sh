#! /usr/bin/env bash
in=/mnt/picea/storage/reference/Populus-tremula/v2.2
out=/mnt/picea/storage/reference/Populus-tremula/v2.2/indices/snpEff
mail=nicolas.delhomme@umu.se
account=facility

# load modules
module load bioinfo-tools snpEff

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

sbatch -A $account --mail-user $mail -o $out/build.out -e $out/build.err \
$(realpath ../../../UPSCb-common/pipeline/runSnpEffBuild.sh) \
$out/snpEff.config Potra02 \
$in/fasta/Potra02_genome.fasta \
$in/gff/Potra02_genes.gff3 \
$in/fasta/Potra02_proteins.fasta \
$in/fasta/Potra02_transcripts.fasta
