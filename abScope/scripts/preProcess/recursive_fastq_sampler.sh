#!/bin/bash

# Verificar se seqkt está instalado
if ! command -v seqtk &> /dev/null; then
    echo "Error: seqtk not found"
    exit 1
fi

# Verifica os argumentos
if [[ $# -lt 2 ]]; then
    echo "Use: $0 <diretorio_input> <diretorio_output> [N=1000]"
    exit 1
fi

INPUT_DIR=$1
OUTPUT_DIR=$2
N=${3:-1000} # Se não fornecido, usa 1000 por padrão

if [-d $OUTPUT_DIR]; then
    echo "Directory $OUTPUT_DIR alredy created, using it."
else 
    # Se não existe, cria
    mkdir -p "$OUTPUT_DIR"
    echo "Directory $OUTPUT_DIR created!"
fi

# Diretório de saída
find "$INPUT_DIR" -type d | while read -r DIR; do
    NEW_DIR="${OUTPUT_DIR}${DIR#$INPUT_DIR}"
    mkdir -p "$NEW_DIR"
done

# Processa os arquivos .fastq
find "$INPUT_DIR" -type f -name "*.fastq" | while read -r FILE; do
    REL_PATH="${FILE#$INPUT_DIR/}"
    OUT_FILE="$OUTPUT_DIR/$REL_PATH"

    echo "Processing: $FILE -> $OUT_FILE"

    # Sampling with seqtk
    seqtk sample "$FILE" "$N" > "$OUT_FILE"
done

echo "Sampling Completed!"
