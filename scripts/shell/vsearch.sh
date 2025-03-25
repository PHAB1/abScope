run_vsearch() {
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
}