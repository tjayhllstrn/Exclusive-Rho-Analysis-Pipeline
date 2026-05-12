#!/bin/bash

set -euo pipefail

DIR="$1"
MAX_JOBS=6   # adjust for your machine
MODEL_PATH="src/gbt/models/model_rga_pass2_inbending"

running=0
mkdir -p logs
source ~/.Exclusive-Rho-Analysis-Pipeline/bin/activate

for file in "$DIR"/nSidis_*.root; do
    echo python3 src/gbt/predict.py "$file" "$MODEL_PATH" "EventTree"
    log="logs/$(basename "$file" .root).log"
    nohup python3 src/gbt/predict.py "$file" "$MODEL_PATH" "EventTree" \
	  > "$log" 2>&1 &

    ((running+=1))
    if (( running >= MAX_JOBS )); then
        wait -n
        ((running-=1))
    fi
done

wait
echo "ALL jobs finished"
