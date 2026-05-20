#!/bin/bash
#SBATCH --partition=production
#SBATCH --account=clas12
#SBATCH --mem=5000
#SBATCH --time=01:00:00
#SBATCH --output=/farm_out/%u/%x-%j-%N.out
#SBATCH --error=/farm_out/%u/%x-%j-%N.err
#SBATCH --job-name=pipeline1_hadd

OUTDIR=$1
FNAME=$2

cd "$OUTDIR"

# Remove existing output file if it exists
[ -f "$FNAME" ] && rm "$FNAME"

# Combine all .root files in OUTDIR into one merged file
hadd -f "$FNAME" ./*.root
