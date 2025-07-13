import pandas as pd
import numpy as np
import glob

# Load configuration
configfile: "config.yaml"

# Output Directories
OUTPUT_DIR="results"
IGCICLE_DIR="results/igcicle"
CLUSTERING_DIR="results/clustering"
CLUSTERING_RAW=f"{CLUSTERING_DIR}/raw"
CLUSTERING_CDR=f"{CLUSTERING_DIR}/cdr"
CONSENSUS_DIR=f"results/consensus"
ANNOTATION_DIR=f"results/annotation"

# log Directory
LOG_DIR="logs"

# Define script and Conda environment based on sequencing type
pre_process_illumina_script = "scripts/preProcess/preProcess_illumina.bash"
pre_process_nano_script = "scripts/preProcess/preProcess_nanopore.bash"
igcicle_script = "scripts/preProcess/igCicle.py"
clusterCdr3_script = "scripts/preProcess/cluster_cdr3.bash"
reGroupCdrCluster_script = "scripts/preProcess/reGroupCdrCluster.py"
map_clusters_script = "scripts/preProcess/mapClusters.py"
polish_script = "scripts/preProcess/polish.py"
define_clones_script = "scripts/script_analysis/recursive_define_clones.sh"
dupcount_script = "scripts/script_analysis/dupcount_col.sh"
list_fasta_script = "scripts/preProcess/list_fasta.py"

# envs
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
        expand(f"{OUTPUT_DIR}/quality_trimmed/{{sample}}/qc.fasta", sample=samples), # quality trimmed
        expand(f"{IGCICLE_DIR}/{{sample}}/VH_seqs.fasta", sample=samples), # IGCICLE VH
        expand(f"{IGCICLE_DIR}/{{sample}}/VL_seqs.fasta", sample=samples),  # IGCICLE VL
        expand(f"{IGCICLE_DIR}/{{sample}}/1_igCicle.tsv", sample=samples), # IGCICLE 1
        expand(f"{IGCICLE_DIR}/{{sample}}/2_igCicle.tsv", sample=samples),  # IGCICLE 2
        expand(f"{CLUSTERING_RAW}/{{sample}}/VH_clusters.uc", sample=samples), # VH clustering
        expand(f"{CLUSTERING_RAW}/{{sample}}/VL_clusters.uc", sample=samples),  # VL clustering
        expand(f"{CLUSTERING_CDR}/{{sample}}/cdr3_clusters/done.txt", sample=samples), # cdr3 mapping
        expand(f"{CLUSTERING_CDR}/{{sample}}/new_clusters/new_cdr3_clusters/done.txt", sample=samples), # cdr3 clustering
        expand(f"{CLUSTERING_CDR}/{{sample}}/new_clusters/new_cdr3_clusters/new_final_clusters", sample=samples), # pre Medaka 
        expand(f"{CLUSTERING_CDR}/{{sample}}/new_clusters/new_cdr3_clusters/new_final_clusters/done.txt", sample=samples), # pre Medaka done file
        expand(f"{CLUSTERING_CDR}/{{sample}}/medaka_consensus", sample=samples),
        expand(f"{CONSENSUS_DIR}/{{sample}}/cdr3_consensus.fasta", sample=samples),
        expand(f"{ANNOTATION_DIR}/{{sample}}/igblastn_annotation.tsv", sample=samples),
        expand(f"{ANNOTATION_DIR}/{{sample}}/igblastn_annotation_dup.tsv", sample=samples)

rule qualityTrimming:
    conda:
        pre_process_env
    input:
        file_1=lambda wildcards: sample_info[wildcards.sample]["file_1"]
    output:
        fq=f"{OUTPUT_DIR}/quality_trimmed/{{sample}}/qc.fastq",
        fa=f"{OUTPUT_DIR}/quality_trimmed/{{sample}}/qc.fasta"
    params:
        file_2=lambda wildcards: sample_info[wildcards.sample]["file_2"],
        generation=lambda wildcards: sample_info[wildcards.sample]["generation"]
    log:
        f"{LOG_DIR}/quality_trimmed/{{sample}}/trimming.log"
    shell:
        """
        # second generation 
        if [[ {params.generation} == "second" ]]; then
            {pre_process_illumina_script} {input.file_1} {params.file_2} {output.fq} >> {log} 2>&1

        # third generation
        else
            NanoFilt {input.file_1} -q 8 -l 700 > {output.fq}
        fi

        # rename the fasta file
        TEMP_FASTA=$(mktemp)
        RENAMED_FASTA="{output.fa}"
        awk 'NR%4==1 {{printf(">%s\\n", substr($0, 2)); next}} NR%4==2 {{print}}' {output.fq} > "$TEMP_FASTA"

        if [[ {config[rename_reads]} == True ]]; then
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
        f"results/quality_trimmed/{{sample}}/qc.fasta"
    output:
        vh=f"{IGCICLE_DIR}/{{sample}}/VH_seqs.fasta",
        vl=f"{IGCICLE_DIR}/{{sample}}/VL_seqs.fasta",
        igcicle1=f"{IGCICLE_DIR}/{{sample}}/1_igCicle.tsv",
        igcicle2=f"{IGCICLE_DIR}/{{sample}}/2_igCicle.tsv",
        cicle1_passed=f"{IGCICLE_DIR}/{{sample}}/cicle_1_passed.fasta",
        cicle2_passed=f"{IGCICLE_DIR}/{{sample}}/cicle_2_passed.fasta"
    params:
        generation=lambda wildcards: sample_info[wildcards.sample]["generation"]
    log:
        f"{LOG_DIR}/igcicle/{{sample}}/igCicle.log"
    shell:
        """
        python3 {igcicle_script} {input} {output.vh} {output.vl} >> {log} 2>&1

	# correct annotations 
	# After this correction:
	# ig_cicle1 and ig_cicle2 corrected #
	ig_db_path="references/database/human/ig_db"

	# Cicle 1
	igblastn -germline_db_J $ig_db_path/IGHLKJ_edit.fasta -germline_db_V $ig_db_path/IGHLKV_edit.fasta -germline_db_D $ig_db_path/IGHD_edit.fasta \
	    -domain_system imgt -query {output.cicle1_passed} -auxiliary_data $ig_db_path/human_gl.aux \
	    -out {output.igcicle1} -num_threads 2 -outfmt 19 -show_translation

	    
	# Cicle 2 if exist
	if [ -f {output.cicle2_passed} ] && [ -s {output.cicle2_passed} ]; then
	    igblastn -germline_db_J $ig_db_path/IGHLKJ_edit.fasta -germline_db_V $ig_db_path/IGHLKV_edit.fasta -germline_db_D $ig_db_path/IGHD_edit.fasta \
	    -domain_system imgt -query {output.cicle2_passed} -auxiliary_data $ig_db_path/human_gl.aux \
	    -out {output.igcicle2} -num_threads 2 -outfmt 19 -show_translation

	else
	    echo "cicle 2 doesn't produced annotation"
	fi
        """

rule rawClustering:
    conda:
        medaka_env
    input:
        vh=f"{IGCICLE_DIR}/{{sample}}/VH_seqs.fasta",
        vl=f"{IGCICLE_DIR}/{{sample}}/VL_seqs.fasta"
    output:
        vh_uc=f"{CLUSTERING_RAW}/{{sample}}/VH_clusters.uc",
        vl_uc=f"{CLUSTERING_RAW}/{{sample}}/VL_clusters.uc",
        vh_cons=f"{CLUSTERING_RAW}/{{sample}}/VH_consensus.fasta",
        vl_cons=f"{CLUSTERING_RAW}/{{sample}}/VL_consensus.fasta"
    log:
        f"{LOG_DIR}/clustering/raw/{{sample}}/clustering.log"
    shell:
        """
	n_lines_vh=$(wc -l < {input.vh})
	if [ -f {input.vh} ] && [ -s {input.vh} ] && [ $n_lines_vh -gt 20 ]; then
            vsearch --cluster_fast {input.vh} --id 0.78 --uc {output.vh_uc} \
                --target_cov 0.9 --minqt 0.9 --consout {output.vh_cons} >> {log} 2>&1
	else
	    touch {output.vh_uc} {output.vh_cons}
	fi
        
	n_lines_vl=$(wc -l < {input.vl})
	if [ -f {input.vl} ] && [ -s {input.vl} ] && [ $n_lines_vl -gt 20 ]; then
            vsearch --cluster_fast {input.vl} --id 0.78 --uc {output.vl_uc} \
                --target_cov 0.9 --minqt 0.9 --consout {output.vl_cons} >> {log} 2>&1
	else
            touch {output.vl_uc} {output.vl_cons}
	fi
        """

rule mapClusters:
    conda:
        vsearch_env
    input:
        vh_fa=f"{CLUSTERING_RAW}/{{sample}}/VH_consensus.fasta",
        vl_fa=f"{CLUSTERING_RAW}/{{sample}}/VL_consensus.fasta",
        vh_uc=f"{CLUSTERING_RAW}/{{sample}}/VH_clusters.uc",
        vl_uc=f"{CLUSTERING_RAW}/{{sample}}/VL_clusters.uc",
        igcicle_1=f"{IGCICLE_DIR}/{{sample}}/1_igCicle.tsv",
        igcicle_2=f"{IGCICLE_DIR}/{{sample}}/2_igCicle.tsv"
    output:
        done=f"{CLUSTERING_CDR}/{{sample}}/cdr3_clusters/done.txt",
        json=f"{CLUSTERING_CDR}/{{sample}}/cdr3_clusters/cluster_seq_id.json", 
        igall=f"{IGCICLE_DIR}/{{sample}}/igAll.csv",
        consensus=directory(f"{CLUSTERING_CDR}/{{sample}}/clusters_consensus")
    log:
        f"{LOG_DIR}/clustering/cdr/{{sample}}/cdr3_clusters/clustering.log"
    shell:
        """
        # process VH ## Arrumar geração do meta_VH e meta_VL aqui ##
	n_lines_vh=$(wc -l < {input.vh_fa})
	if [ -f {input.vh_fa} ] && [ -s {input.vh_fa} ] && [ $n_lines_vh -gt 20 ]; then # depois substituir pelo n minimo de seqs p/ formar um cluster
            python3 {map_clusters_script} {input.vh_fa} {input.vh_uc} {input.igcicle_1} \
                {input.igcicle_2} VH {output.done} >> {log} 2>&1 
	else
	    echo "empty {wildcards.sample} for VH"
	fi

        # process VL
	n_lines_vl=$(wc -l < {input.vl_fa})
	if [ -f {input.vl_fa} ] && [ -s {input.vl_fa} ] && [ $n_lines_vl -gt 20 ]; then
            python3 {map_clusters_script} {input.vl_fa} {input.vl_uc} {input.igcicle_1} \
                {input.igcicle_2} VL {output.done} >> {log} 2>&1
        else
	    echo "empty {wildcards.sample} for VL"
	fi

        output_dir=$(dirname {output.done})
	python3 {list_fasta_script} $output_dir {output.done}
        """

rule cdrClustering:
    conda:
        vsearch_env
    input:
        file_list=f"{CLUSTERING_CDR}/{{sample}}/cdr3_clusters/done.txt",
    output:
        done=f"{CLUSTERING_CDR}/{{sample}}/new_clusters/new_cdr3_clusters/done.txt",
    log:
        f"{LOG_DIR}/clustering/cdr/{{sample}}/new_clusters/new_cdr3_clusters/cdr_clustering.log"
    shell:
        """
        {clusterCdr3_script} {input.file_list} {output.done} >> {log} 2>&1
        """

# prepare medaka files (clusters) to polish with Medaka
rule prePolish:
    conda:
        vsearch_env
    input:
        json=f"{CLUSTERING_CDR}/{{sample}}/cdr3_clusters/cluster_seq_id.json",
        new_clust_list=f"{CLUSTERING_CDR}/{{sample}}/new_clusters/new_cdr3_clusters/done.txt",
        igall=f"{IGCICLE_DIR}/{{sample}}/igAll.csv"
    output:
        nfc=directory(f"{CLUSTERING_CDR}/{{sample}}/new_clusters/new_cdr3_clusters/new_final_clusters"),
        done=f"{CLUSTERING_CDR}/{{sample}}/new_clusters/new_cdr3_clusters/new_final_clusters/done.txt",
        json_out=f"{CLUSTERING_CDR}/{{sample}}/new_clusters/new_cdr3_clusters/new_final_clusters/cluster_seq_id.json"
    log:
        f"{LOG_DIR}/clustering/cdr/{{sample}}/new_clusters/new_cdr3_clusters/new_final_clusters/cd3_clusterization.log"
    shell:
        """
        # VH type 
        python3 {reGroupCdrCluster_script} -j {input.json} -l {input.new_clust_list} -i {input.igall} -c VH -nfc {output.nfc} -jo {output.json_out} >> {log} 2>&1
         
        # VL type 
        python3 {reGroupCdrCluster_script} -j {input.json} -l {input.new_clust_list} -i {input.igall} -c VL -nfc {output.nfc} -jo {output.json_out} >> {log} 2>&1

        # create the done file
        output_dir=$(dirname {output.done})
        #ls $output_dir/*.fasta >> {output.done}
	python3 {list_fasta_script} $output_dir {output.done}
        """

# run medaka on selected groups
rule Polish:
    conda:
        medaka_env
    input:
        CLUSTERS_CONSENSUS_DIR=f"{CLUSTERING_CDR}/{{sample}}/clusters_consensus", 
        NEWS_FINAL_CLUSTERS=f"{CLUSTERING_CDR}/{{sample}}/new_clusters/new_cdr3_clusters/new_final_clusters/done.txt",
        json=f"{CLUSTERING_CDR}/{{sample}}/new_clusters/new_cdr3_clusters/new_final_clusters/cluster_seq_id.json"
    output:
        MEDAKA_OUTPUT_DIR=directory(f"{CLUSTERING_CDR}/{{sample}}/medaka_consensus"),
        FINAL_CONSENSUS=f"{CONSENSUS_DIR}/{{sample}}/cdr3_consensus.fasta"
    shell:
        """
	python3 {polish_script} {input.CLUSTERS_CONSENSUS_DIR} {input.NEWS_FINAL_CLUSTERS} {input.json} 2 {output.FINAL_CONSENSUS} {output.MEDAKA_OUTPUT_DIR}
	"""

# Medaka Final cicle --| colocar depois 

# Annotate the consensus sequences with igblast
rule igAnnotateConsensus:
    conda:
        igblast_env
    input:
        consensus=f"{CONSENSUS_DIR}/{{sample}}/cdr3_consensus.fasta",
    output:
        annot=f"{ANNOTATION_DIR}/{{sample}}/igblastn_annotation.tsv",
        clones=f"{ANNOTATION_DIR}/{{sample}}/igblastn_annotation_clone-pass.tsv",
        dup_count=f"{ANNOTATION_DIR}/{{sample}}/igblastn_annotation_dup.tsv"
    log:
        f"{LOG_DIR}/annotation/{{sample}}/annotation.log"
    shell:
        """
        igblastn -germline_db_J references/database/human/ig_db/IGHLKJ_edit.fasta -germline_db_V references/database/human/ig_db/IGHLKV_edit.fasta -germline_db_D references/database/human/ig_db/IGHD_edit.fasta \
         -query {input.consensus} -outfmt 19 \
         -show_translation -auxiliary_data references/database/human/ig_db/human_gl.aux \
         -num_threads 8 -out {output.annot} \
         -domain_system imgt
        
	# define clones script
	output_dir=$(dirname {output.annot})
	bash {define_clones_script} $output_dir $output_dir

	# add duplicate count column 
        {dupcount_script} {output.clones} {output.dup_count}
        """ 
