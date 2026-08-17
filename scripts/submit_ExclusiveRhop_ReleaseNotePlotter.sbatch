#!/bin/bash

#SBATCH --partition=production
#SBATCH --account=clas12
#SBATCH --mem-per-cpu=2000
#SBATCH --nodes=1
#SBATCH --chdir=/scratch/slurm
#SBATCH --output=/farm_out/%u/%x-%j-%N.out
#SBATCH --error=/farm_out/%u/%x-%j-%N.err
#SBATCH --job-name=ExclusiveRhop_ReleaseNotePlotter
#SBATCH --time=3:00:00

cd /w/hallb-scshelf2102/clas12/users/tjhellst/Exclusive-Rho-Analysis-Pipeline/
source setup.sh

root -l -b -q 'scripts/ExclusiveRhop_ReleaseNotePlotter.cpp("/volatile/clas12/users/tjhellst/ExclusiveRhoPlus_RGA_processedData/pippi0_merged_in_pass2/pippi0_merged_in_pass2.root",true)'

