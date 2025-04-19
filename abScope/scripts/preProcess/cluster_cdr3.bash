#!/bin/bash

# Verifica se os argumentos foram passados corretamente
if [[ $# -ne 2 ]]; then
    echo "Uso: $0 <input_paths> <output_done>"
    exit 1
fi

paths="$1"
output="$2"
output_dir=$(dirname "${output}")

while IFS= read -r path; do
    filename=$(basename "${path}" ".fasta")
    vsearch --cluster_fast "${path}" --id 0.85 --target_cov 0.95 --minseqlength 6 --threads "${threads}" \
        --uc "${output_dir}/${filename}.uc" --consout "${output_dir}/${filename}.fasta"
    ls "${output_dir}/${filename}.fasta" >> "${output}"
done < "${paths}"