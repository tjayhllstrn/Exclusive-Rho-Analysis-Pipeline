#!/bin/bash

#SBATCH --partition=production
#SBATCH --account=clas12
#SBATCH --mem-per-cpu=5000
#SBATCH --nodes=2
#SBATCH --array=0-119
#SBATCH --chdir=/scratch/slurm
#SBATCH --output=/farm_out/%u/%x-%j-%N.out
#SBATCH --error=/farm_out/%u/%x-%j-%N.err
#SBATCH --job-name=pippi0_builder_spring2019
#SBATCH --time=3:00:00

pwd
cd /w/hallb-scshelf2102/clas12/users/tjhellst/Exclusive-Rho-Analysis-Pipeline/
pwd
source /u/home/tjhellst/.Exclusive-Rho-Analysis-Pipeline/bin/activate

DIR=$1
INPUT_FILE_BASE=nSidis_*.root

# Enable nullglob to get empty array if no matches (instead of literal pattern)
shopt -s nullglob
FILES=($DIR/$INPUT_FILE_BASE) 
shopt -u nullglob

# Debug output
echo "INPUT_DIR: $DIR"
echo "INPUT_FILE_BASE: $INPUT_FILE_BASE"
echo "Pattern: $DIR/$INPUT_FILE_BASE"
echo "Number of files found: ${#FILES[@]}"
echo "Task ID: $SLURM_ARRAY_TASK_ID"

# Check if files were found
if [ ${#FILES[@]} -eq 0 ]; then
    echo "ERROR: No files matched pattern $DIR/$INPUT_FILE_BASE"
    exit 1
fi


echo "Running job ${SLURM_ARRAY_TASK_ID[@]} on file: ${FILES[$SLURM_ARRAY_TASK_ID]}"

clas12root -l -b -q "macros/pippi0Builder.C(\"${FILES[$SLURM_ARRAY_TASK_ID]}\")"
