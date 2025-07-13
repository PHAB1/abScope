#!/bin/bash

input="$1"
output="$2"

# Opcional: Adicionar verificação de argumentos
if [ -z "$input" ] || [ -z "$output" ]; then
  echo "Uso: $0 <arquivo_de_entrada> <arquivo_de_saida>"
  exit 1
fi

echo "Arquivo de entrada: $input"
echo "Arquivo de saída: $output"

awk 'BEGIN{FS="\t"; OFS="\t"}
NR==1 {print $0, "duplicate_count"; next} # Adiciona o header na primeira linha e pula para a próxima
{
  dupcount_value = 1
  if (match($0, /DUPCOUNT=([0-9]+)/, arr)) {
    dupcount_value = arr[1]
  }
  print $0, dupcount_value
}' "$input" > "$output"

echo "Processamento concluído."
