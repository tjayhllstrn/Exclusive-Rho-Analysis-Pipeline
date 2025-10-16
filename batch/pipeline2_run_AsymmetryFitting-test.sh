#!/bin/bash

#SBATCH --partition=production
#SBATCH --account=clas12
#SBATCH --mem-per-cpu=3000
#SBATCH --nodes=2
#SBATCH --chdir=/scratch/slurm
#SBATCH --output=/farm_out/%u/%x-%j-%N.out
#SBATCH --error=/farm_out/%u/%x-%j-%N.err
#SBATCH --job-name=pipeline2_AsymFit
#SBATCH --time=1:00:00
#./pipeline2_run_AsymmetryFitting-test.sh pippi0_RGAinbending_zbinning /w/hallb-scshelf2102/clas12/users/tjhellst/clas-ana-scaffold-tyler/out/test/nSidis_005032.root /w/hallb-scshelf2102/clas12/users/tjhellst/clas-ana-scaffold-tyler/out/test/nSidis_005036.root
CONFIG=$1
shift
INPUT_FILES=("$@")

cd /w/hallb-scshelf2102/clas12/users/tjhellst/clas-ana-scaffold-tyler

./macros/AsymmetryFitting.py "$CONFIG" "${INPUT_FILES[@]}"
