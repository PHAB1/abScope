import pandas as pd

# Load configuration
configfile: "config.yaml"

# Directories
FASTA_DIR = "quality_trimmed"
IGCICLE_DIR = "results/igCicle"
CLUSTERING_DIR = "results/clustering"
LOG_DIR = "logs"

# Coleta todos os arquivos .fasta no diretório
FASTA_FILES = glob_wildcards(os.path.join(FASTA_DIR, "{sample}.fasta")).sample

# Define script and Conda environment based on sequencing type
pre_process_illumina_script = "scripts/preProcess/preProcess_illumina.bash"
pre_process_nano_script = "scripts/preProcess/preProcess_nanopore.bash"
pre_process_env = "envs/preProcess.yaml"
igblast_env = "envs/igblast.yaml"
vsearch_env = "envs/vsearch.yaml"

rule all:
    input:
        "results/quality_trimmed",
        IGCICLE_DIR,
        CLUSTERING_DIR

rule quality_trimming:
    conda:
        pre_process_env
    input:
        samples=config["input"]
    output:
        directory("results/quality_trimmed")
    log:
        f"{LOG_DIR}/log_trimming.txt"
    shell:
        """
        # illumina pre process
        bash {pre_process_illumina_script} {input.samples} {output} >> {log} 2>&1

        # nanopore pre process
        bash {pre_process_nano_script} {input.samples} {output} >> {log} 2>&1
        """

rule igAnnotate:
    conda:
        igblast_env
    input:
        fastas="results/quality_trimmed"
    output:
        directory(IGCICLE_DIR)
    log:
        f"{LOG_DIR}/log_igBlast.txt"
    shell:
        """
        for fasta in {input.fastas}/*.fasta; do
            python scripts/preProcess/igCicle.py "$fasta" {output} >> {log} 2>&1
        done
        """

rule rawClustering:
    conda:
        vsearch_env
    input:
        samples=IGCICLE_DIR
    output:
        directory(CLUSTERING_DIR)
    log:
        f"{LOG_DIR}/log_clustering.txt"
    shell:
        """
        run_vsearch() {{
            local input_file="$1"
            local output_dir="$2"
            local chain="$3"

            vsearch \
                --cluster_fast "$input_file" \
                --id 0.75 \
                --uc "$output_dir/${{chain}}_clusters.uc" \
                --target_cov 0.9 \
                --minqt 0.9 \
                --consout "$output_dir/${{chain}}_consensus.fasta" 
        }}

        for dir in {input.samples}/*/; do
            basedir=$(basename "$dir")


            # Executando para VH e VL
            mkdir -p "{output}/${{basedir}}"
            run_vsearch "${{dir}}VH_seqs.fasta" "{output}/${{basedir}}" VH
            run_vsearch "${{dir}}VL_seqs.fasta" "{output}/${{basedir}}" VL
            
        done
        """