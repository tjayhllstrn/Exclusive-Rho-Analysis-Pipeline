#!/bin/bash

#SBATCH --partition=production
#SBATCH --account=clas12
#SBATCH --mem-per-cpu=3000
#SBATCH --nodes=2
#SBATCH --array=0-0 
#SBATCH --chdir=/scratch/slurm
#SBATCH --output=/farm_out/%u/%x-%j-%N.out
#SBATCH --error=/farm_out/%u/%x-%j-%N.err
#SBATCH --job-name=pipeline0_runSingleHipo_pippi0
#SBATCH --time=5:00:00

cd /w/hallb-scshelf2102/clas12/users/tjhellst/Exclusive-Rho-Analysis-Pipeline/
source /u/home/tjhellst/.Exclusive-Rho-Analysis-Pipeline/bin/activate
module load clas12root

INPUT_DIR=$1
OUTPUT_DIR=$2
FILES=($INPUT_DIR/nSidis_*.hipo)

TYPE=pippi0

echo "Running job ${SLURM_ARRAY_TASK_ID[@]} on file: ${FILES[$SLURM_ARRAY_TASK_ID]}"

./scripts/run_single_hipo.rb --type $TYPE --input ${FILES[$SLURM_ARRAY_TASK_ID]} --outdir $OUTPUT_DIR -n -1

ˇ