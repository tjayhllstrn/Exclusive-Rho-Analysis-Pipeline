#!/bin/bash

set -euo pipefail

IN_DIR="$1"
ASYM_TYPE="$2"
ASYM_TYPE_BKG="$3"
MAX_JOBS=6   # adjust for your machine

running=0

mkdir -p logs

for file in "$IN_DIR"/clasdis*.root; do
    echo clas12root -l -b -q "macros/InjectionStudy.cpp(\"$file\", \"$ASYM_TYPE\", \"$ASYM_TYPE_BKG\")"
    log="logs/$(basename "$file" .root).log"
    nohup clas12root -l -b -q "macros/InjectionStudy.cpp(\"$file\", \"$ASYM_TYPE\", \"$ASYM_TYPE_BKG\")" \
	  > "$log" 2>&1 &

    ((running+=1))
    if (( running >= MAX_JOBS )); then
        wait -n
        ((running-=1))
    fi
done

wait
echo "ALL jobs finished"
