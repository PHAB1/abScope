import pandas as pd
import numpy as np
import glob

# Load configuration
configfile: "config.yaml"

# Directories
<<<<<<< HEAD
OUTPUT_DIR="results"
IGCICLE_DIR="results/igCicle"
CLUSTERING_DIR="results/clustering"
LOG_DIR="logs"
=======
FASTA_DIR = "quality_trimmed"
IGCICLE_DIR = "results/igCicle"
CLUSTERING_DIR = "results/clustering"
LOG_DIR = "logs"

# Coleta todos os arquivos .fasta no diretório
FASTA_FILES = glob_wildcards(os.path.join(FASTA_DIR, "{sample}.fasta")).sample
>>>>>>> 154e883a95c4fb7a3db782670acfcde579c6ab53

# Define script and Conda environment based on sequencing type
pre_process_illumina_script = "scripts/preProcess/preProcess_illumina.bash"
pre_process_nano_script = "scripts/preProcess/preProcess_nanopore.bash"
pre_process_env = "envs/preProcess.yaml"
igblast_env = "envs/igblast.yaml"
vsearch_env = "envs/vsearch.yaml"
<<<<<<< HEAD
cdr_clustering_env = "envs/cdrClusters.yaml"

# get the list of columns from the input file
SAMPLES_DF=pd.read_csv(config["input"], sep=",")

# Create a dictionary to map sample names to their corresponding columns
sample_info=SAMPLES_DF.set_index('samples').T.to_dict()
samples=SAMPLES_DF['samples'].tolist()
ch_types=SAMPLES_DF['ch_type'].tolist()

rule all:
    input:
        [f"{OUTPUT_DIR}/quality_trimmed/{sample}/{ch_type}/teste.fasta" for sample, ch_type in zip(samples, ch_types)]
        #IGCICLE_DIR
        #expand(f"{CLUSTERING_DIR}/{{sample}}/VH_clusters.uc", sample=samples),
        #expand(f"{CLUSTERING_DIR}/{{sample}}/VL_clusters.uc", sample=samples)
=======

rule all:
    input:
        "results/quality_trimmed",
        IGCICLE_DIR,
        CLUSTERING_DIR
>>>>>>> 154e883a95c4fb7a3db782670acfcde579c6ab53

rule quality_trimming:
    conda:
        pre_process_env
    input:
        lambda wildcards: sample_info[wildcards.sample]["file"]
    output:
        f"{OUTPUT_DIR}/quality_trimmed/{{sample}}/{{ch_type}}/teste.fasta"
    shell:
        """
        # second generation ## Agora é necessário passar as variáveis novamente, arrume os arquivos correspondentes --##
        bash scripts/preProcess/preProcess_illumina.bash \ 
            {input} {output} >> {log} 2>&1

        # third generation
        bash scripts/preProcess/preProcess_nanopore.bash \
            {input} {output} >> {log} 2>&1        
        """

# inativo em um momento    
rule igAnnotate:
    conda:
        igblast_env
    input:
        fastas="results/quality_trimmed"
    output:
<<<<<<< HEAD
=======
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
>>>>>>> 154e883a95c4fb7a3db782670acfcde579c6ab53
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
<<<<<<< HEAD
        vh=f"{IGCICLE_DIR}/{{sample}}/VH_seqs.fasta",
        vl=f"{IGCICLE_DIR}/{{sample}}/VL_seqs.fasta"
    output:
        vh_uc=f"{CLUSTERING_DIR}/{{sample}}/VH_clusters.uc",
        vl_uc=f"{CLUSTERING_DIR}/{{sample}}/VL_clusters.uc",
        vh_cons=f"{CLUSTERING_DIR}/{{sample}}/VH_consensus.fasta",
        vl_cons=f"{CLUSTERING_DIR}/{{sample}}/VL_consensus.fasta"
    log:
        f"{LOG_DIR}/clustering_{{sample}}.txt"
    shell:
        """
        vsearch --cluster_fast {input.vh} --id 0.75 --uc {output.vh_uc} \
            --target_cov 0.9 --minqt 0.9 --consout {output.vh_cons}
        
        vsearch --cluster_fast {input.vl} --id 0.75 --uc {output.vl_uc} \
            --target_cov 0.9 --minqt 0.9 --consout {output.vl_cons}
=======
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
>>>>>>> 154e883a95c4fb7a3db782670acfcde579c6ab53
        """