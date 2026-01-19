#!/bin/bash

#SBATCH --partition=production
#SBATCH --account=clas12
#SBATCH --mem-per-cpu=4000
#SBATCH --nodes=4
#SBATCH --chdir=/scratch/slurm
#SBATCH --output=/farm_out/%u/%x-%j-%N.out
#SBATCH --error=/farm_out/%u/%x-%j-%N.err
#SBATCH --job-name=pipeline3_AsymFit
#SBATCH --time=2:00:00

FILE_NAME=$1
FILE_DIR=$2
FIT_TYPE=$3
CONFIG=$4

cd /w/hallb-scshelf2102/clas12/users/tjhellst/Exclusive-Rho-Analysis-Pipeline/

clas12root -l -b -q "macros/AsymmetryFitting.cpp(\"$FILE_NAME\", \"$FILE_DIR\", \"$FIT_TYPE\", \"$CONFIG\")"
