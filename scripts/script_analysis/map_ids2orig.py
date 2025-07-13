import csv
import os
import sys

def substituir_ids(arquivo_mapeamento, arquivo_alvo):
    """
    Substitui ocorrências de IDs em um arquivo alvo usando um mapeamento de IDs.

    Args:
        arquivo_mapeamento (str): Caminho para o arquivo CSV de mapeamento (orig, new).
        arquivo_alvo (str): Caminho para o arquivo a ter os IDs substituídos.
    """
    mapeamento = {}
    try:
        with open(arquivo_mapeamento, 'r', newline='') as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                if len(row) == 2:
                    orig, new = row
                    mapeamento[new] = orig
                else:
                    print(f"Aviso: Linha mal formatada no arquivo de mapeamento: {row}")
    except FileNotFoundError:
        print(f"Erro: Arquivo de mapeamento '{arquivo_mapeamento}' não encontrado.")
        sys.exit(1)

    nome_base_alvo = os.path.splitext(os.path.basename(arquivo_alvo))[0]
    extensao_alvo = os.path.splitext(arquivo_alvo)[1]
    arquivo_saida = f"{nome_base_alvo}_substituido{extensao_alvo}"

    try:
        with open(arquivo_alvo, 'r') as infile, open(arquivo_saida, 'w') as outfile:
            for linha in infile:
                nova_linha = linha
                for novo_id, original_id in mapeamento.items():
                    nova_linha = nova_linha.replace(novo_id, original_id)
                outfile.write(nova_linha)
    except FileNotFoundError:
        print(f"Erro: Arquivo alvo '{arquivo_alvo}' não encontrado.")
        sys.exit(1)

    print(f"Substituições completas. O novo arquivo é: {arquivo_saida}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Uso: python script.py <arquivo_mapeamento.csv> <arquivo_alvo>")
        sys.exit(1)
    arquivo_mapeamento = sys.argv[1]
    arquivo_alvo = sys.argv[2]
    substituir_ids(arquivo_mapeamento, arquivo_alvo)
