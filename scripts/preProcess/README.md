# vsearch for clustering 
vsearch --cluster_fast ../../nanopore/minIONLinear_Data/barcode10_VL_quality-pass.fasta --id 0.9 --centroids quality-passbarcode10_VL_quality_vsearch.fasta --uc clusters.uc --target_cov 0.9 --minqt 0.9 --consout quality-passbarcode10_VL_quality_consensus.fasta

# verify if sequences in each cluster has at least 98% identity (to do, use igblastn for this task)


# verify Snps frequencies between each cluster and split cluster if frequence > value [val. to discover] (to do, use bcftools for this task)


# Best consensus with medaka


# presto + igblastn rescursive pipelines, with dup_counts and cprimer=None included


## nanopore processing -> Old bash

#!/bin/bash

# vsearch for clustering 
#vsearch --cluster_fast ../../nanopore/minIONLinear_Data/barcode06_VL_quality-pass.fasta --id 0.9 --centroids quality-passbarcode06_VL_quality_vsearch.fasta --uc clusters.uc --target_cov 0.9 --minqt 0.9 --consout quality-passbarcode06_VL_quality_consensus.fasta

# verify if sequences in each cluster has at least 98% identity (to do, use igblastn for this task)


# verify Snps frequencies between each cluster and split cluster if frequence > value [val. to discover] (to do, use bcftools for this task)


# Best consensus with medaka


# presto + igblastn rescursive pipelines, with dup_counts and cprimer=None included


## New steps ##
# use samtools fastq input.bam > output.fastq ; to transform bam to fastq

# Quality
NanoFilt barcode06.fastq -q 8 -l 2800 > barcode06_testeFiltered.fastq

# Fastq to Fasta
awk 'NR%4==1 {printf(">%s\n", substr($0, 2)); next} NR%4==2 {print}' barcode06_testeFiltered.fastq > barcode06_Filtered.fasta

# Rename sequences ids 
python scripts/fasta_rename_id.py barcode06_Filtered.fasta barcode06_Filtered_renamed.fasta

# igCicle for recursively retain ig regions
python scripts/igCicle.py barcode06_Filtered_renamed.fasta

# VH and VL igCicle files
# VH
igblastn -germline_db_J ig_db/IGHLKJ_edit.fasta -germline_db_V ig_db/IGHLKV_edit.fasta -germline_db_D ig_db/IGHD_edit.fasta -query igCicle/VH_seqs.fasta -outfmt 19 -show_translation -auxiliary_data ig_db/human_gl.aux -num_threads 8 -out igCicle/VH_seqs.tsv -domain_system kabat

# VL
igblastn -germline_db_J ig_db/IGHLKJ_edit.fasta -germline_db_V ig_db/IGHLKV_edit.fasta -germline_db_D ig_db/IGHD_edit.fasta -query igCicle/VL_seqs.fasta -outfmt 19 -show_translation -auxiliary_data ig_db/human_gl.aux -num_threads 8 -out igCicle/VL_seqs.tsv -domain_system kabat

# vsearch VH
#vsearch --cluster_fast ../../nanopore/minIONLinear_Data/barcode06_VL_quality-pass.fasta --id 0.9 --centroids quality-passbarcode06_VL_quality_vsearch.fasta --uc clusters.uc --target_cov 0.9 --minqt 0.9 --consout quality-passbarcode06_VL_quality_consensus.fasta
#vsearch --cluster_fast igCicle/VH_seqs.fasta --id 0.9 --centroids barcode14_vsearch.fasta --uc clusters.uc --target_cov 0.9 --minqt 0.9 --consout barcode06_vsearch_consensus.fasta
vsearch --cluster_size igCicle/VH_seqs.fasta --id 0.75 --uc clusters.uc --target_cov 0.9 --minqt 0.9 --consout barcode06_vsearch_consensus.fasta

# cluster by CDR3
python scripts/createCDR3Clusters.py "barcode06" "VH"

# medaka consensus using consensus as draft for cicle 2
#medaka_consensus -i igCicle/VH_seqs.fasta -d barcode06_processed.fasta -o medaka_cicle2_VH -f -x
medaka_consensus -i igCicle/VH_seqs.fasta -d barcode06_processed.fasta -o medaka_cicle2_VH -f -x -g -t 8

# igblastn on final consensus file
igblastn -germline_db_J ig_db/IGHLKJ_edit.fasta -germline_db_V ig_db/IGHLKV_edit.fasta -germline_db_D ig_db/IGHD_edit.fasta -query barcode06_processed.fasta -outfmt 19 -show_translation -auxiliary_data ig_db/human_gl.aux -num_threads 8 -out barcode06_procesed_VH.tsv -domain_system kabat

# vsearch VL
vsearch --cluster_size igCicle/VL_seqs.fasta --id 0.75 --uc clusters.uc --target_cov 0.9 --minqt 0.9 --consout barcode06_vsearch_consensus.fasta

# cluster by CDR3
python scripts/createCDR3Clusters.py "barcode06" "VL"

# medaka consensus using consensus as draft for cicle 2
#medaka_consensus -i igCicle/VL_seqs.fasta -d barcode06_processed.fasta -o medaka_cicle2_VL -f -x
medaka_consensus -i igCicle/VL_seqs.fasta -d barcode06_processed.fasta -o medaka_cicle2_VL -f -x -g -t 8

# igblastn on final consensus file
igblastn -germline_db_J ig_db/IGHLKJ_edit.fasta -germline_db_V ig_db/IGHLKV_edit.fasta -germline_db_D ig_db/IGHD_edit.fasta -query barcode06_processed.fasta -outfmt 19 -show_translation -auxiliary_data ig_db/human_gl.aux -num_threads 8 -out barcode06_procesed_VL.tsv -domain_system kabat
