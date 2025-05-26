import pandas as pd
import numpy as np # Importar numpy para usar np.nan

def calculate_region_identity(query_segment, germline_segment):
    """
    Calcula a porcentagem de identidade entre dois segmentos de sequência alinhados.
    A identidade é calculada como (Número de Matches) / (Comprimento do Segmento) * 100.
    Matches são posições onde os caracteres são idênticos.
    O comprimento do segmento é o número total de posições no segmento alinhado.
    """
    if len(query_segment) != len(germline_segment):
        # Isso não deveria acontecer com segmentos extraídos de alinhamentos de mesmo comprimento
        print(f"Aviso: segmentos com comprimentos diferentes ({len(query_segment)} vs {len(germline_segment)})")
        return 0.0 # Ou np.nan, dependendo de como quer tratar isso

    total_length = len(query_segment)
    if total_length == 0:
        return 0.0

    # Contar matches (posições onde os caracteres são idênticos)
    matches = sum(q == g for q, g in zip(query_segment, germline_segment))

    identity = (matches / total_length) * 100

    return identity

def add_region_identities(input_filepath, output_filepath):
    """
    Lê o output tabular do igblastn, calcula a identidade para as regiões CDR e FR,
    e salva a tabela atualizada.
    """
    try:
        # Ler o arquivo de entrada. Assumimos que é separado por tabulação.
        # low_memory=False é útil para arquivos grandes para evitar warnings sobre tipos mistos.
        df = pd.read_csv(input_filepath, sep='\t', low_memory=False)
    except FileNotFoundError:
        print(f"Erro: Arquivo de entrada não encontrado em {input_filepath}")
        return
    except Exception as e:
        print(f"Erro ao ler o arquivo de entrada: {e}")
        return

    # Definir as regiões e seus nomes de colunas correspondentes para início/fim e identidade.
    # Use os nomes das colunas exatamente como aparecem no seu arquivo.
    regions = {
        'fwr1': ('fwr1_start', 'fwr1_end', 'fwr1_identity'),
        'cdr1': ('cdr1_start', 'cdr1_end', 'cdr1_identity'),
        'fwr2': ('fwr2_start', 'fwr2_end', 'fwr2_identity'),
        'cdr2': ('cdr2_start', 'cdr2_end', 'cdr2_identity'),
        'fwr3': ('fwr3_start', 'fwr3_end', 'fwr3_identity'),
        'cdr3': ('cdr3_start', 'cdr3_end', 'cdr3_identity'), # Note que CDR3 é geralmente calculado de forma diferente, mas aqui usamos as coordenadas fornecidas.
        'fwr4': ('fwr4_start', 'fwr4_end', 'fwr4_identity'),
    }

    # Verificar se as colunas essenciais existem no DataFrame
    required_cols = ['sequence_alignment', 'germline_alignment']
    for _, (start_col, end_col, _) in regions.items():
        required_cols.append(start_col)
        required_cols.append(end_col)

    if not all(col in df.columns for col in required_cols):
        missing = [col for col in required_cols if col not in df.columns]
        print(f"Erro: O arquivo de entrada não possui as colunas necessárias: {', '.join(missing)}")
        return

    # Calcular a identidade para cada região
    # Usamos itertuples para iteração eficiente sobre as linhas
    new_data = []
    for row in df.itertuples(index=True): # index=True inclui o índice como o primeiro elemento da tupla
        row_data = row._asdict() # Converte a tupla de volta para um dicionário para fácil acesso pelas colunas

        for region_name, (start_col, end_col, identity_col) in regions.items():
            row_data[identity_col] = 0.0 # Inicializa a nova coluna de identidade para esta linha

            try:
                # Obter os valores de início e fim. Usamos getattr para acessar colunas por nome de string.
                start = getattr(row, start_col)
                end = getattr(row, end_col)
                seq_align = getattr(row, 'sequence_alignment')
                germ_align = getattr(row, 'germline_alignment')

                # Verificar se as coordenadas são válidas e os alinhamentos não são vazios/NaN
                if pd.notna(start) and pd.notna(end) and start > 0 and end >= start and pd.notna(seq_align) and pd.notna(germ_align):
                    # Ajustar as coordenadas para fatiamento (slicing) baseado em 0 no Python
                    # As coordenadas do igblastn são geralmente baseadas em 1 e inclusivas.
                    start_idx = int(start) - 1
                    end_idx = int(end) # O slicing no Python é exclusivo no final, então o fim inclusivo vira o índice exclusivo

                    query_segment = seq_align[start_idx:end_idx]
                    germline_segment = germ_align[start_idx:end_idx]

                    # Verificar se os segmentos têm o mesmo comprimento antes de calcular a identidade
                    if len(query_segment) == len(germline_segment):
                        identity = calculate_region_identity(query_segment, germline_segment)
                        row_data[identity_col] = identity
                    else:
                         # Isso pode indicar um problema com as coordenadas ou os alinhamentos
                         print(f"Aviso: Diferença no comprimento do segmento para a linha de índice {row.Index}, região {region_name}. Pulando o cálculo para esta região nesta linha.")
                         row_data[identity_col] = np.nan # Define como NaN se houver problema

                else:
                    # Coordenadas inválidas ou faltando, definir identidade como 0.0 ou NaN
                    row_data[identity_col] = 0.0 # Ou np.nan se preferir indicar dado faltante

            except Exception as e:
                print(f"Erro ao processar a linha de índice {row.Index}, região {region_name}: {e}")
                row_data[identity_col] = np.nan # Em caso de erro inesperado, define como NaN

        new_data.append(row_data)

    # Criar um novo DataFrame a partir dos dados processados
    updated_df = pd.DataFrame(new_data)

    # Garantir que a ordem das colunas originais seja mantida, adicionando as novas colunas no final
    original_columns = df.columns.tolist()
    new_identity_columns = [regions[region][2] for region in regions]
    updated_df = updated_df[original_columns + new_identity_columns]


    # Salvar o DataFrame atualizado em um novo arquivo
    try:
        # Usar sep='\t' para salvar como arquivo separado por tabulação
        # index=False para não escrever o índice do DataFrame como uma coluna
        updated_df.to_csv(output_filepath, sep='\t', index=False)
        print(f"Identidades das regiões adicionadas com sucesso. Output salvo em {output_filepath}")
    except Exception as e:
        print(f"Erro ao salvar o arquivo de output: {e}")

# --- Configuração ---
# ALtere o nome do arquivo de entrada para o seu arquivo .tsv ou .txt
input_file = "all_phage_airr_c_call_primer.tsv"
# Altere o nome do arquivo de saída para o nome desejado para o novo arquivo
output_file = "identidades_de_regiao.tsv"
# ---------------------

if __name__ == "__main__":
    add_region_identities(input_file, output_file)
