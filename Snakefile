import pandas as pd
import numpy as np
import glob

# Load configuration
configfile: "config.yaml"

# Directories
OUTPUT_DIR="results"
IGCICLE_DIR="results/igCicle"
CLUSTERING_DIR="results/clustering"
LOG_DIR="logs"

# Define script and Conda environment based on sequencing type
pre_process_illumina_script = "scripts/preProcess/preProcess_illumina.bash"
pre_process_nano_script = "scripts/preProcess/preProcess_nanopore.bash"
pre_process_env = "envs/preProcess.yaml"
igblast_env = "envs/igblast.yaml"
vsearch_env = "envs/vsearch.yaml"
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
        """