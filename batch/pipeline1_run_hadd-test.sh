#!/bin/bash
#SBATCH --partition=production
#SBATCH --account=clas12
#SBATCH --mem=2500
#SBATCH --time=00:30:00
#SBATCH --output=/farm_out/%u/%x-%j-%N.out
#SBATCH --error=/farm_out/%u/%x-%j-%N.err
#SBATCH --job-name=pipeline1_hadd
#pipeline1_run_hadd-test.sh /w/hallb-scshelf2102/clas12/users/tjhellst/clas-ana-scaffold-tyler/out/test hadd_test.root

OUTDIR=$1
FNAME=$2

cd "$OUTDIR"

# Combine all .root files in OUTDIR into one merged file
hadd -f "$FNAME" ./*.root