import os
import uuid
import argparse
from Bio import SeqIO
import pandas as pd

def rename_fasta_ids(input_fasta, output_fasta):
    # Abrir o arquivo de entrada
    meta_file = "%s.meta"%("".join(output_fasta.split(".")[0]))
    meta_dict = {
            "orig":[],
            "new":[]
            }
    with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
        # Iterar sobre cada sequência no arquivo FASTA
        for record in SeqIO.parse(infile, "fasta"):
            # Gerar um novo ID aleatório
            new_id = str(uuid.uuid4())
            # record meta id
            meta_dict["orig"].append(record.id)
            meta_dict["new"].append(new_id)
            # Atribuir o novo ID ao record
            record.id = new_id
            record.description = new_id  # Remove a descrição original
            # Escrever a sequência no novo arquivo FASTA
            outfile.write(f">{record.id}\n{str(record.seq)}\n")
            
    meta_df = pd.DataFrame(meta_dict)
    meta_df.to_csv(meta_file, index=False)

if __name__ == "__main__":
    # Configurando o argparse para aceitar argumentos de linha de comando
    parser = argparse.ArgumentParser(description="Renomear IDs de um arquivo FASTA para IDs aleatórios.")
    parser.add_argument("input_fasta", help="Caminho para o arquivo FASTA de entrada.")
    parser.add_argument("output_fasta", help="Caminho para o arquivo FASTA de saída com IDs renomeados.")
    
    args = parser.parse_args()

    # Chamar a função com os argumentos fornecidos
    rename_fasta_ids(args.input_fasta, args.output_fasta)
