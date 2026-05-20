#!/bin/bash

#SBATCH --partition=production
#SBATCH --account=clas12
#SBATCH --mem-per-cpu=5000
#SBATCH --nodes=2
#SBATCH --array=0-352
#SBATCH --chdir=/scratch/slurm
#SBATCH --output=/farm_out/%u/%x-%j-%N.out
#SBATCH --error=/farm_out/%u/%x-%j-%N.err
#SBATCH --job-name=Inj_Asym_RGAMCpass2
#SBATCH --time=0:30:00

cd /w/hallb-scshelf2102/clas12/users/tjhellst/Exclusive-Rho-Analysis-Pipeline/

INDIR=$1
INJASYM_TYPE=$2
BKGASYM_TYPE=$3
INPUT_FILE_BASE=clasdis_*.root

# Enable nullglob to get empty array if no matches (instead of literal pattern)
shopt -s nullglob
FILES=($INDIR/$INPUT_FILE_BASE) 
shopt -u nullglob

# Debug output
echo "INPUT_DIR: $INDIR"
echo "INPUT_FILE_BASE: $INPUT_FILE_BASE"
echo "Pattern: $INDIR/$INPUT_FILE_BASE"
echo "Number of files found: ${#FILES[@]}"
echo "Task ID: $SLURM_ARRAY_TASK_ID"

# Check if files were found
if [ ${#FILES[@]} -eq 0 ]; then
    echo "ERROR: No files matched pattern $INDIR/$INPUT_FILE_BASE"
    exit 1
fi


echo "Running job ${SLURM_ARRAY_TASK_ID[@]} on file: ${FILES[$SLURM_ARRAY_TASK_ID]}"

clas12root -l -b -q "macros/InjectionStudy.cpp(\"${FILES[$SLURM_ARRAY_TASK_ID]}\", \"$INJASYM_TYPE\", \"$BKGASYM_TYPE\")"
