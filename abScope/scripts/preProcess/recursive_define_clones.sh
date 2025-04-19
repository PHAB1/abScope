#!/bin/bash

# Access input directory
echo "input directory: $1"
echo ""

# Acess Output directory
echo "Output directory: $2"
echo ""

for file in $1/*.tsv; do
    filename=$(basename "$file")
    filename_no_ext=${filename%.*}
    echo "Processing: $filename_no_ext"
    DefineClones.py -d "$file" --outdir "$2" --format airr --nproc 12 --mode gene --act first
done
