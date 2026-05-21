#!/bin/bash

#SBATCH --partition=production
#SBATCH --account=clas12
#SBATCH --mem-per-cpu=5000
#SBATCH --nodes=2
#SBATCH --chdir=/scratch/slurm
#SBATCH --output=/farm_out/%u/%x-%j-%N.out
#SBATCH --error=/farm_out/%u/%x-%j-%N.err
#SBATCH --job-name=MC_data_loss_allMCfa18_inb
#SBATCH --time=5:00:00

pwd
cd /w/hallb-scshelf2102/clas12/users/tjhellst/Exclusive-Rho-Analysis-Pipeline/
pwd

MC_file=$1
Data_file=$2


clas12root -l -b -q "macros/MC_data_loss.cpp(\"$MC_file\",\"$Data_file\")"