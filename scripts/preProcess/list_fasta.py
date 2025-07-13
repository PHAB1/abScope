import os
import sys

def listar_fastas(diretorio_saida, arquivo_output):
    """
    Lista todos os arquivos .fasta em um diretório e salva seus caminhos em um arquivo.

    Args:
        diretorio_saida (str): O caminho do diretório a ser pesquisado.
        arquivo_output (str): O caminho do arquivo onde a lista de FASTA será salva.
    """
    try:
        with open(arquivo_output, 'w') as f_out:
            for item in os.listdir(diretorio_saida):
                # Constrói o caminho completo para o item
                caminho_completo = os.path.join(diretorio_saida, item)
                # Verifica se é um arquivo e termina com .fasta
                if os.path.isfile(caminho_completo) and item.endswith(".fasta"):
                    f_out.write(caminho_completo + '\n')
        print(f"Arquivos .fasta listados e salvos em: {arquivo_output}")
    except OSError as e:
        sys.stderr.write(f"Erro ao acessar o diretório ou escrever no arquivo: {e}\n")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write("Uso: python seu_script.py <diretorio_saida> <arquivo_output.done>\n")
        sys.exit(1)

    output_dir = sys.argv[1]
    output_done_file = sys.argv[2]
    listar_fastas(output_dir, output_done_file)
