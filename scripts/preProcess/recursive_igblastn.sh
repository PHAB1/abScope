#!/bin/bash

# Access input directory
echo "input directory: $1"
echo ""

# Acess Output directory
echo "Output directory: $2"
echo ""

# Check if at least one argument is provided
if [ $# -ge 3 ]; then
    prefix="$3"
else
    prefix=""
fi

ig_db_path="ig_db/"

for file in $1/*.fasta; do
    filename=$(basename "$file")
    filename_no_ext=${filename%.*}
    echo "Processing: $filename_no_ext"
    echo "Output file: $prefix$filename_no_ext.tsv"
    igblastn -germline_db_J "$ig_db_path/IGHLKJ_edit.fasta" -germline_db_V "$ig_db_path/IGHLKV_edit.fasta" -germline_db_D "$ig_db_path/IGHD_edit.fasta" -query "$file" -outfmt 19 -show_translation -auxiliary_data "$ig_db_path"/human_gl.aux -num_threads 12 -out "$2/$prefix$filename_no_ext.tsv"
    #igblastn -germline_db_J ../ig_db/IGHLKJ_edit.fasta -germline_db_V ../ig_db/IGHLKV_edit.fasta -germline_db_D ../ig_db/IGHD_edit.fasta -query "$file" -outfmt 19 -show_translation -auxiliary_data ../ig_db/human_gl.aux -num_threads 12 -out "$2/$prefix$filename_no_ext.tsv"
done
