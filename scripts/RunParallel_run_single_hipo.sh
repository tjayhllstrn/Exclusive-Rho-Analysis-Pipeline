#!/bin/bash

set -euo pipefail

INDIR="$1"
OUTDIR="$2"
MAX_JOBS=6   # adjust for your machine

running=0

mkdir -p logs

source ~/.Exclusive-Rho-Analysis-Pipeline/bin/activate

for file in "$INDIR"/nSidis_0050*.hipo; do
    echo ruby scripts/run_single_hipo.rb --type pippi0 --input "$file" --outdir "$OUTDIR" -n -1
    log="logs/$(basename "$file" .hipo).log"
    nohup ruby scripts/run_single_hipo.rb --type pippi0 --input "$file" --outdir "$OUTDIR" -n -1 \
	  > "$log" 2>&1 &

    ((running+=1))
    if (( running >= MAX_JOBS )); then
        wait -n
        ((running-=1))
    fi
done

wait
echo -e "\e[1;32mALL jobs finished\e[0m"
