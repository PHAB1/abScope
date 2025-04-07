#!/bin/bash

# Verifica se os argumentos foram passados corretamente
if [[ $# -ne 2 ]]; then
    echo "Uso: $0 <csv_samples_path> <output_dir>"
    exit 1
fi

# Define as variáveis a partir dos argumentos
CSV_FILE="$1"
OUTPUT_DIR="$2"

# Verifica se o arquivo CSV existe
if [[ ! -f "$CSV_FILE" ]]; then
    echo "Erro: Arquivo CSV '$CSV_FILE' não encontrado!"
    exit 1
fi

# Cria o diretório de saída se não existir
mkdir -p "$OUTPUT_DIR"

# Lê o CSV e processa cada arquivo FASTQ
awk -F, 'NR>1 {print $1,$5}' "$CSV_FILE" | while read -r INPUT_FASTQ generation; do
    
    if [ "$generation" != "third" ]; then
        continue
    fi

    # Verifica se o arquivo FASTQ existe
    if [[ ! -f "$INPUT_FASTQ" ]]; then
        echo "Erro: Arquivo FASTQ '$INPUT_FASTQ' não encontrado! Pulando..."
        continue
    fi

    # Define nomes de saída baseados no nome do arquivo de entrada
    BASENAME=$(basename "$INPUT_FASTQ" .fastq)
    FILTERED_FASTQ="$OUTPUT_DIR/${BASENAME}_filtered.fastq"
    FILTERED_FASTA="$OUTPUT_DIR/${BASENAME}_filtered.fasta"
    RENAMED_FASTA="$OUTPUT_DIR/${BASENAME}_filtered_renamed.fasta"

    echo "Processando: $INPUT_FASTQ"

    # Filtragem de qualidade
    NanoFilt "$INPUT_FASTQ" -q 8 -l 2800 > "$FILTERED_FASTQ"

    # Conversão para FASTA
    TEMP_FASTA=$(mktemp)
    awk 'NR%4==1 {printf(">%s\n", substr($0, 2)); next} NR%4==2 {print}' "$FILTERED_FASTQ" > "$TEMP_FASTA"

    # Renomear sequências
    python scripts/preProcess/fasta_rename_id.py "$TEMP_FASTA" "$RENAMED_FASTA"
    rm "$TEMP_FASTA"
done
