#!/bin/bash
#SBATCH --partition=production
#SBATCH --account=clas12
#SBATCH --mem=3500
#SBATCH --time=4:00:00
#SBATCH --output=/farm_out/%u/%x-%j-%N.out
#SBATCH --error=/farm_out/%u/%x-%j-%N.err
#SBATCH --job-name=pipeline2_hadd_merge

# Get all arguments except the last one (input files)
# Last argument is the output directory
ARGS=("$@")
OUTDIR="${ARGS[-1]}"
unset 'ARGS[-1]'
INPUT_FILES=("${ARGS[@]}")

# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"

# Output filename: OUTDIR/basename(OUTDIR).root
OUTDIR_BASENAME=$(basename "$OUTDIR")
OUTPUT_FILE="$OUTDIR/${OUTDIR_BASENAME}.root"

echo "Merging ${#INPUT_FILES[@]} files into $OUTPUT_FILE"
echo "Input files: ${INPUT_FILES[@]}"

# Combine all input .root files into one merged file
hadd -f "$OUTPUT_FILE" "${INPUT_FILES[@]}"