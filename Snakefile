import pandas as pd

# Load configuration
configfile: "config.yaml"

# Validate sequencing type
if config["sequencing"] not in ["illumina", "nanopore"]:
    raise ValueError("Error: 'sequencing' must be either 'illumina' or 'nanopore'.")

# Define script and Conda environment based on sequencing type
if config["sequencing"] == "illumina":
    pre_process_script = "scripts/preProcess/preProcess_illumina.bash"
    pre_process_env = "envs/fastp.yaml"
elif config["sequencing"] == "nanopore":
    pre_process_script = "scripts/preProcess/preProcess_nanopore.bash"
    pre_process_env = "envs/nanopore_preProcess.yaml"

rule all:
    input:
        "quality_trimmed"

rule quality_trimming:
    conda:
        pre_process_env
    input:
        samples=config["input"]
    output:
        directory("quality_trimmed")
    log:
        "logs/log_trimming.txt"
    shell:
        """
        bash {pre_process_script} {input.samples} {output} >> {log} 2>&1
        """

