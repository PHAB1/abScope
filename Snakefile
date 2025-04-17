import pandas as pd
import numpy as np
import glob

# Load configuration
configfile: "config.yaml"

# Directories
OUTPUT_DIR="results"
IGCICLE_DIR="results/igcicle"
CLUSTERING_DIR="results/clustering"
CLUSTERING_RAW=f"{CLUSTERING_DIR}/raw"
CLUSTERING_CDR=f"{CLUSTERING_DIR}/cdr"
LOG_DIR="logs"

# Define script and Conda environment based on sequencing type
pre_process_illumina_script = "scripts/preProcess/preProcess_illumina.bash"
pre_process_nano_script = "scripts/preProcess/preProcess_nanopore.bash"
igcicle_script = "scripts/preProcess/igCicle.py"
#cdrCluster_script = "scripts/preProcess/createCDR3Clusters.py"
map_clusters = "scripts/preProcess/mapClusters.py"
pre_process_env = "envs/preProcess.yaml"
igblast_env = "envs/igblast.yaml"
vsearch_env = "envs/vsearch.yaml"
medaka_env = "envs/medaka.yaml"

# get the list of columns from the input file
SAMPLES_DF=pd.read_csv(config["input"], sep=",")

# Create a dictionary to map sample names to their corresponding columns
sample_info=SAMPLES_DF.set_index('samples').T.to_dict()
samples=SAMPLES_DF['samples'].tolist()
ch_types=SAMPLES_DF['ch_type'].tolist()

rule all:
    input:
        [f"{OUTPUT_DIR}/quality_trimmed/{sample}/{ch_type}/qc.fasta" for sample, ch_type in zip(samples, ch_types)],
        [f"{IGCICLE_DIR}/{sample}/{ch_type}/VH_seqs.fasta" for sample, ch_type in zip(samples, ch_types)], # IGCICLE VH
        [f"{IGCICLE_DIR}/{sample}/{ch_type}/VL_seqs.fasta" for sample, ch_type in zip(samples, ch_types)],  # IGCICLE VL
        [f"{IGCICLE_DIR}/{sample}/{ch_type}/1_igCicle.tsv" for sample, ch_type in zip(samples, ch_types)], # IGCICLE 1
        [f"{IGCICLE_DIR}/{sample}/{ch_type}/2_igCicle.tsv" for sample, ch_type in zip(samples, ch_types)],  # IGCICLE 2
        [f"{CLUSTERING_RAW}/{sample}/{ch_type}/VH_clusters.uc" for sample, ch_type in zip(samples, ch_types)], # VH clustering
        [f"{CLUSTERING_RAW}/{sample}/{ch_type}/VL_clusters.uc" for sample, ch_type in zip(samples, ch_types)],  # VL clustering
        [f"{CLUSTERING_CDR}/{sample}/{ch_type}/cdr3_clusters/done_vh.txt" for sample, ch_type in zip(samples, ch_types)], # cdr3 mapping
        [f"{CLUSTERING_CDR}/{sample}/{ch_type}/cdr3_clusters/done_vl.txt" for sample, ch_type in zip(samples, ch_types)], # cdr3 mapping
        [f"{CLUSTERING_CDR}/{sample}/{ch_type}/new_clusters/new_cdr3_clusters/done_vh.txt" )], # cdr3 VH clustering
        [f"{CLUSTERING_CDR}/{sample}/{ch_type}/new_clusters/new_cdr3_clusters/done_vl.txt" )] # cdr3 VL clustering
        #[f"{CLUSTERING_CDR}/{sample}/{ch_type}/VH_final_consensus.fasta" for sample, ch_type in zip(samples, ch_types)], # cdr3 re-clustering + consensus
        #[f"{CLUSTERING_CDR}/{sample}/{ch_type}/VL_final_consensus.fasta" for sample, ch_type in zip(samples, ch_types)] # cdr3 re-clustering + consensus

rule qualityTrimming:
    conda:
        pre_process_env
    input:
        file_1=lambda wildcards: sample_info[wildcards.sample]["file_1"]
    output:
        fq=f"{OUTPUT_DIR}/quality_trimmed/{{sample}}/{{ch_type}}/qc.fastq",
        fa=f"{OUTPUT_DIR}/quality_trimmed/{{sample}}/{{ch_type}}/qc.fasta"
    params:
        file_2=lambda wildcards: sample_info[wildcards.sample]["file_2"],
        generation=lambda wildcards: sample_info[wildcards.sample]["generation"]
    log:
        f"{LOG_DIR}/quality_trimmed/{{sample}}/{{ch_type}}/trimming.log"
    shell:
        """
        # second generation 
        if [[ {params.generation} == "second" ]]; then
            {pre_process_illumina_script} {input.file_1} {params.file_2} {output.fq} >> {log} 2>&1

        # third generation
        else
            NanoFilt {input.file_1} -q 8 -l 2800 > {output.fq}
        fi

        # rename the fasta file
        TEMP_FASTA=$(mktemp)
        RENAMED_FASTA="{output.fa}"
        awk 'NR%4==1 {{printf(">%s\\n", substr($0, 2)); next}} NR%4==2 {{print}}' {output.fq} > "$TEMP_FASTA"

        if [[ {config[rename_reads]} == "True" ]]; then
            python scripts/preProcess/fasta_rename_id.py "$TEMP_FASTA" "{output.fa}"
        else
            mv "$TEMP_FASTA" "{output.fa}"
        fi
        """

# run IgCicle (multiple igblastn cicles)
rule igAnnotate:
    conda:
        igblast_env
    input:
        f"results/quality_trimmed/{{sample}}/{{ch_type}}/qc.fasta"
    output:
        vh=f"{IGCICLE_DIR}/{{sample}}/{{ch_type}}/VH_seqs.fasta",
        vl=f"{IGCICLE_DIR}/{{sample}}/{{ch_type}}/VL_seqs.fasta",
        igcicle1=f"{IGCICLE_DIR}/{{sample}}/{{ch_type}}/1_igCicle.tsv",
        igcicle2=f"{IGCICLE_DIR}/{{sample}}/{{ch_type}}/2_igCicle.tsv"
    params:
        generation=lambda wildcards: sample_info[wildcards.sample]["generation"]
    log:
        f"{LOG_DIR}/igcicle/{{sample}}/{{ch_type}}/igCicle.log"
    shell:
        """
        python3 {igcicle_script} {input} {output.vh} {output.vl} >> {log} 2>&1
        """

rule rawClustering:
    conda:
        medaka_env
    input:
        vh=f"{IGCICLE_DIR}/{{sample}}/{{ch_type}}/VH_seqs.fasta",
        vl=f"{IGCICLE_DIR}/{{sample}}/{{ch_type}}/VL_seqs.fasta"
    output:
        vh_uc=f"{CLUSTERING_RAW}/{{sample}}/{{ch_type}}/VH_clusters.uc",
        vl_uc=f"{CLUSTERING_RAW}/{{sample}}/{{ch_type}}/VL_clusters.uc",
        vh_cons=f"{CLUSTERING_RAW}/{{sample}}/{{ch_type}}/VH_consensus.fasta",
        vl_cons=f"{CLUSTERING_RAW}/{{sample}}/{{ch_type}}/VL_consensus.fasta"
    log:
        f"{LOG_DIR}/clustering/raw/{{sample}}/{{ch_type}}/clustering.log"
    shell:
        """
        vsearch --cluster_fast {input.vh} --id 0.75 --uc {output.vh_uc} \
            --target_cov 0.9 --minqt 0.9 --consout {output.vh_cons} >> {log} 2>&1
        
        vsearch --cluster_fast {input.vl} --id 0.75 --uc {output.vl_uc} \
            --target_cov 0.9 --minqt 0.9 --consout {output.vl_cons} >> {log} 2>&1
        """

rule mapClusters:
    conda:
        medaka_env
    input:
        vh_fa=f"{CLUSTERING_RAW}/{{sample}}/{{ch_type}}/VH_consensus.fasta",
        vl_fa=f"{CLUSTERING_RAW}/{{sample}}/{{ch_type}}/VL_consensus.fasta",
        vh_uc=f"{CLUSTERING_RAW}/{{sample}}/{{ch_type}}/VH_clusters.uc",
        vl_uc=f"{CLUSTERING_RAW}/{{sample}}/{{ch_type}}/VL_clusters.uc",
        igcicle_1=f"{IGCICLE_DIR}/{{sample}}/{{ch_type}}/1_igCicle.tsv",
        igcicle_2=f"{IGCICLE_DIR}/{{sample}}/{{ch_type}}/2_igCicle.tsv"
    output:
        vh=f"{CLUSTERING_CDR}/{{sample}}/{{ch_type}}/cdr3_clusters/done_vh.txt",
        vl=f"{CLUSTERING_CDR}/{{sample}}/{{ch_type}}/cdr3_clusters/done_vl.txt"
    log:
        f"{LOG_DIR}/clustering/cdr/{{sample}}/{{ch_type}}/clustering.log"
    shell:
        """
        # process VH
        python3 {map_clusters} {input.vh_fa} {input.vh_uc} {input.igcicle_1} \
            {input.igcicle_2} VH {output.vh} >> {log} 2>&1 

        # process VL
        python3 {map_clusters} {input.vl_fa} {input.vl_uc} {input.igcicle_1} \
            {input.igcicle_2} VL {output.vl} >> {log} 2>&1

        touch {output.vh}
        touch {output.vl}
        """

rule cdrClustering:
    conda:
        vsearch_env
    input:
        vh=f"{CLUSTERING_CDR}/{{sample}}/{{ch_type}}/cdr3_clusters/{id1}.fasta",
        vl=f"{CLUSTERING_CDR}/{{sample}}/{{ch_type}}/cdr3_clusters/{id2}.fasta"
    output:
        vh=f"{CLUSTERING_CDR}/{{sample}}/{{ch_type}}/new_clusters/new_cdr3_clusters/done_vh.txt"
        vh=f"{CLUSTERING_CDR}/{{sample}}/{{ch_type}}new_clusters/new_cdr3_clusters/done_vl.txt"
    shell:
        """
        # process VH
        vsearch --cluster_fast {input.vh} --id 0.85 --target_cov 0.9 --minseqlength 6 --uc %s.uc --consout *.fasta

        # process VL
        vsearch --cluster_fast {input.vl} --id 0.85 --target_cov 0.9 --minseqlength 6 --uc %s.uc --consout *.fasta
	"""
        

#rule cdrClustering:
#    conda:
#        medaka_env
#    input:
#        vh_fa=f"{CLUSTERING_RAW}/{{sample}}/{{ch_type}}/VH_consensus.fasta",
#        vl_fa=f"{CLUSTERING_RAW}/{{sample}}/{{ch_type}}/VL_consensus.fasta",
#       vh_uc=f"{CLUSTERING_RAW}/{{sample}}/{{ch_type}}/VH_clusters.uc",
#       vl_uc=f"{CLUSTERING_RAW}/{{sample}}/{{ch_type}}/VL_clusters.uc",
#        igcicle_1=f"{IGCICLE_DIR}/{{sample}}/{{ch_type}}/1_igCicle.tsv",
#        igcicle_2=f"{IGCICLE_DIR}/{{sample}}/{{ch_type}}/2_igCicle.tsv"
#    output:
#        vh=f"{CLUSTERING_CDR}/{{sample}}/{{ch_type}}/VH_final_consensus.fasta",
#        vl=f"{CLUSTERING_CDR}/{{sample}}/{{ch_type}}/VL_final_consensus.fasta"
#    log:
#        f"{LOG_DIR}/clustering/cdr/{{sample}}/{{ch_type}}/clustering.log"
#    shell:
#        """
#        # process VH
#        python3 {cdrCluster_script} {input.vh_fa} {input.vh_uc} {input.igcicle_1} \
#            {input.igcicle_2} VH {output.vh} 
#
#        # process VL
#        python3 {cdrCluster_script} {input.vl_fa} {input.vl_uc} {input.igcicle_1} \
#            {input.igcicle_2} VL {output.vl} 
#        """
