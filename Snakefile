import pandas as pd

# Load configuration
configfile: "config.yaml"

# Directories
FASTA_DIR = "quality_trimmed"

# Coleta todos os arquivos .fasta no diretório
FASTA_FILES = glob_wildcards(os.path.join(FASTA_DIR, "{sample}.fasta")).sample

# Define script and Conda environment based on sequencing type
pre_process_illumina_script = "scripts/preProcess/preProcess_illumina.bash"
pre_process_nano_script = "scripts/preProcess/preProcess_nanopore.bash"
pre_process_env = "envs/preProcess.yaml"
igblast_env = "envs/igblast.yaml"

rule all:
    input:
        "results/quality_trimmed",
        "results/igCicle"

rule quality_trimming:
    conda:
        pre_process_env
    input:
        samples=config["input"]
    output:
        directory("results/quality_trimmed")
    log:
        "logs/log_trimming.txt"
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
        fastas = directory("results/quality_trimmed")
    output:
        directory("results/igCicle")
    log:
        "logs/log_igBlast.txt"
    shell:
        """
        for fasta in {input.fastas}/*.fasta; do
            python scripts/preProcess/igCicle.py "$fasta"
        done
        """