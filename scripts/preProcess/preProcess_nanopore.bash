#!/bin/bash

# Verifica se os argumentos foram passados corretamente
if [[ $# -ne 2 ]]; then
    echo "Uso: $0 <input_fastq> <output_dir>"
    exit 1
fi

# Define as variáveis a partir dos argumentos
INPUT_FASTQ="$1"
OUTPUT_DIR="$2"

# Verifica se o arquivo FASTQ existe
if [[ ! -f "$INPUT_FASTQ" ]]; then
    echo "Error: FASTQ file path: '$INPUT_FASTQ' not found!"
    exit 1
fi

# Cria o diretório de saída se não existir
mkdir -p "$OUTPUT_DIR"

# Filtragem de qualidade
FILTERED_FASTQ="$OUTPUT_DIR/filtered.fastq"
NanoFilt "$INPUT_FASTQ" -q 8 -l 2800 > "$FILTERED_FASTQ"

# Conversão para FASTA
FILTERED_FASTA="$OUTPUT_DIR/filtered.fasta"
awk 'NR%4==1 {printf(">%s\n", substr($0, 2)); next} NR%4==2 {print}' "$FILTERED_FASTQ" > "$FILTERED_FASTA"

# Renomear sequências
RENAMED_FASTA="$OUTPUT_DIR/filtered_renamed.fasta"
python scripts/fasta_rename_id.py "$FILTERED_FASTA" "$RENAMED_FASTA"
