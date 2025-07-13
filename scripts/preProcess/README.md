# `abScope`: Analysis Pipeline for Immune Repertoire and Phage Display

`abScope` is a Snakemake-based computational pipeline for processing and analyzing antibody repertoire and phage display library data. It is designed to handle both short-read (e.g., Illumina) and long-read (e.g., Nanopore) sequencing data.

## Core Features
* **Immune Repertoire Processing**: Analyzes diversity, clonal expansion, and VDJ gene usage from B-cell data using an R/Shiny app for visual exploration and filtering of the processed repertoire data.
* **Long-Read Error Correction**: Includes a dedicated module to correct high-error rate sequences from technologies like Oxford Nanopore, producing high-quality consensus sequences.

## Prerequisites
1.  **Snakemake**: To execute the pipeline.
2.  **Conda**: For automatic management of software environments and dependencies.
3.  **R & Shiny**: Required to run the interactive analysis interface.
    ```r
    # In an R console
    install.packages("shiny")
    ```

## How to Use

### 1. Prepare the Input File (`samples.csv`)

The main input is a CSV file that specifies the location and metadata for each sequencing file. It must contain the header: `file_1,file_2,samples,group,generation,ch_type`.

* **`file_1`**: Path to the primary FASTQ file (R1 for paired-end, or the single file for single-end/long-reads).
* **`file_2`**: Path to the R2 FASTQ file for paired-end data. **Leave this column blank for single-end or long-read data.**
* **`samples`**: A unique name/identifier for the sample.
* **`group`**: The experimental group (e.g., `R0`, `R3`, `vaccinated`).
* **`generation`**: Sequencing technology. Use `second` for short-reads or `third` for long-reads.
* **`ch_type`**: if short reads, the antibody chain type being analyzed (e.g., `VH`, `VL`).

**Example `samples.csv`:**
```csv
file_1,file_2,samples,group,generation,ch_type
/path/to/Amostra1_S1_L001_R1_001.fastq,/path/to/Amostra1_S1_L001_R2_001.fastq,Amostra1_illu_VH,R0,second,VH
/path/to/Amostra2_S2_L001_R1_001.fastq,/path/to/Amostra2_S2_L001_R2_001.fastq,Amostra2_illu_VH,R0,second,VH
/path/to/Amostra4_S4_L001_R1_001.fastq,,Amostra4_illu_VL,R0,second,VL
/path/to/nanopore_reads.fastq,,Sample_long_read,R3,third,VH
```

### 2. Configure and Run the Pipeline

**Configuration**: Edit the `config.yaml` file to point to your `samples.csv` and adjust other analysis parameters.

**Execution**: From the project's main directory, run the following command. The pipeline will use Conda to create isolated software environments for each step, ensuring reproducibility.
```bash
snakemake --use-conda --cores <number_of_cores>
```

### 3. Run the Interactive Analysis (Shiny App)

After processing your data, use the Shiny app for interactive analysis of the results.

1.  Open R or RStudio.
2.  Run the following command to launch the application:
    ```r
    shiny::runApp('scripts/immcApp')
    ```
