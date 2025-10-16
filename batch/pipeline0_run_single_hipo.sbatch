#!/bin/bash

#SBATCH --partition=production
#SBATCH --account=clas12
#SBATCH --mem-per-cpu=2500
#SBATCH --nodes=2
#SBATCH --array=0-0 
#SBATCH --chdir=/scratch/slurm
#SBATCH --output=/farm_out/%u/%x-%j-%N.out
#SBATCH --error=/farm_out/%u/%x-%j-%N.err
#SBATCH --job-name=pipeline0_runSingleHipo
#SBATCH --time=0:30:00

cd /w/hallb-scshelf2102/clas12/users/tjhellst/clas-ana-scaffold-tyler

INPUT_DIR=$1
OUTPUT_DIR=$2
FILES=($INPUT_DIR/nSidis_*)

echo "Running job ${SLURM_ARRAY_TASK_ID[@]} on file: ${FILES[$SLURM_ARRAY_TASK_ID]}"

/w/hallb-scshelf2102/clas12/users/tjhellst/clas-ana-scaffold-tyler/scripts/run_single_hipo.rb --type pippi0 
    --input ${FILES[$SLURM_ARRAY_TASK_ID]} 
    --outdir $OUTPUT_DIR -n -1

