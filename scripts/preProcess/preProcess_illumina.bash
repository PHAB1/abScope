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
    echo "Error: CSV file path: '$CSV_FILE' ; not found!"
    exit 1
fi

# Cria o diretório de saída se não existir
mkdir -p "$OUTPUT_DIR"

# Lê o CSV e processa cada amostra
awk -F, 'NR>1 {print $1, $2, $3}' "$CSV_FILE" | while read -r file sample paired; do
    if [[ "$paired" == "1" ]]; then
        r1_file="$file"
    elif [[ "$paired" == "2" ]]; then
        r2_file="$file"
    fi

    # Quando tivermos os dois arquivos (R1 e R2), executamos o fastp
    if [[ -n "$r1_file" && -n "$r2_file" ]]; then
        merged_out="$OUTPUT_DIR/${sample}_merged.fastq"
        no_overlap_r1="$OUTPUT_DIR/${sample}_NoOverlap_R1.fastq"
        no_overlap_r2="$OUTPUT_DIR/${sample}_NoOverlap_R2.fastq"
        r2_lowqual="$OUTPUT_DIR/${sample}_R2_lowQual.fastq"
        r1_lowqual="$OUTPUT_DIR/${sample}_R1_lowQual.fastq"
        html="$OUTPUT_DIR/${sample}.html"
        json="$OUTPUT_DIR/${sample}.json"

        echo "Executando fastp para $sample..."
        fastp -i "$r1_file" -I "$r2_file" --merge --merged_out "$merged_out" \
              -q 15 --out1 "$no_overlap_r1" --out2 "$no_overlap_r2" \
              --unpaired1 "$r2_lowqual" --unpaired2 "$r1_lowqual" \
	      --html "$html" --json "$json"

        # Resetar variáveis para evitar que o próximo par use arquivos errados
        unset r1_file r2_file
    fi
done
