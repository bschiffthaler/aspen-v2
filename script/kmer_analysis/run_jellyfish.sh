#!/bin/bash -l
#SBATCH -n 8
#SBATCH --mem=64G
#SBATCH -t 08:00:00

set -euxo pipefail

JF=/mnt/picea/projects/singularity/jellyfish-2.3.0.simg

# Mount all picea into singularity to prevent io errors
MOUNTBASE='/mnt/picea'

OUTBASE='/mnt/picea/projects/aspseq/nstreet/kmer_ml/sex/count/'

# Gets the longest shared substring
PREFIX=$(sharedPrefix.py "$1" "$2")
PREFIX=$(basename "$PREFIX")

# Doesn't have to be accurate, just consistent
GENOMESIZE=450000000

OUT="${OUTBASE}${PREFIX}_counts.jf"
COUNT="${OUTBASE}${PREFIX}_counts.jf.cov"

if [ ! -f $COUNT ]; then
  cat $1 $2 | pigz -p 8 -d -c | awk '{ if(NR % 4 == 2) { sum += length($0) } } END { print sum }' \
    > $COUNT
fi

# Calculate KMER bounds
N=$(cat $COUNT)

MIN=$(awk -v N=$N -v G=$GENOMESIZE 'BEGIN { printf "%.0f", N / G / 3 }')

# Minimum is always at least 1
if [ $MIN -lt 1 ]; then
  MIN=1
fi

singularity run \
  -B ${MOUNTBASE}:${MOUNTBASE}\
  $JF count -F 2 -s 8G -m 31 -C -L $MIN -o $OUT -t 8 \
    <(gzip -d -c $1 | awk '{if (NR %4 == 1) {print ">"$1} else if (NR %4 ==2) print}') \
    <(gzip -d -c $2 | awk '{if (NR %4 == 1) {print ">"$1} else if (NR %4 ==2) print}')
