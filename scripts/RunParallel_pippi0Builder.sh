#!/bin/bash

set -euo pipefail

DIR="$1"
MAX_JOBS=6   # adjust for your machine

running=0

mkdir -p logs

for file in "$DIR"/*.root; do
    echo clas12root -l -b -q "macros/pippi0Builder.C(\"$file\")"
    log="logs/$(basename "$file" .root).log"
    nohup clas12root -l -b -q "macros/pippi0Builder.C(\"$file\")" \
	  > "$log" 2>&1 &

    ((running+=1))
    if (( running >= MAX_JOBS )); then
        wait -n
        ((running-=1))
    fi
done

wait
echo "ALL jobs finished"
