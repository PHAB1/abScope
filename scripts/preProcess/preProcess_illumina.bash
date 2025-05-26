#!/bin/bash

# Define Inputs
file_1="$1"
file_2="$2"

# define Outputs
OUTPUT="$3"
OUTPUT_DIR=$(dirname "$OUTPUT")
html="$OUTPUT_DIR/fp.html"
json="$OUTPUT_DIR/fp.json"

# Case single end
if [[ -z "$file_2" || "$file_2" == "nan" ]]; then
    lowqual="$OUTPUT_DIR/samplowQual.fastq"

    fastp -i "$file_1" -q 15 -o "$OUTPUT" \
            --html "$html" --json "$json"

# case paired end
else
    merged_out_fasta="$OUTPUT_DIR/${sample}_joined.fasta"
    no_overlap_r1="$OUTPUT_DIR/NoOverlap_R1.fastq"
    no_overlap_r2="$OUTPUT_DIR/NoOverlap_R2.fastq"
    r1_lowqual="$OUTPUT_DIR/R1_lowQual.fastq"
    r2_lowqual="$OUTPUT_DIR/R2_lowQual.fastq"

    fastp -i "$file_1" -I "$file_2" --merge --merged_out "$OUTPUT" \
            -q 15 --out1 "$no_overlap_r1" --out2 "$no_overlap_r2" \
            --unpaired1 "$r2_lowqual" --unpaired2 "$r1_lowqual" \
        --html "$html" --json "$json"
fi
