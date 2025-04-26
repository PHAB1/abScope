import shutil
import subprocess
import glob
import pandas as pd
from Bio import SeqIO
import concurrent.futures
import argparse
import json
import os

parser = argparse.ArgumentParser(
    description="Run Medaka consensus on clustered sequences."
)
parser.add_argument(
    "clusters_consensus_dir",
    help="Path to the directory containing cluster consensus sequences.",
)
parser.add_argument(
    "news_final_clusters",
    help="Path to the directory containing individual sequence files.",
)
parser.add_argument(
    "jsonPath",
    help="Path to json id mapping.",
)
parser.add_argument("threads", help="Threads Number")
parser.add_argument(
    "final_consensus",
    help="Path to the final output fasta files.",
)
parser.add_argument("medaka_output_dir", help="Path to the output directory.")
args = parser.parse_args()

clusters_consensus_dir = args.clusters_consensus_dir
news_final_clusters = args.news_final_clusters 
threads = 8
medaka_output = args.medaka_output_dir
final_consensus = args.final_consensus

output_dir = os.path.dirname(medaka_output)
final_consensus = open(final_consensus, "w")

#cluster_seq_id = pd.DataFrame(cluster_seq_id)
with open(args.jsonPath, 'r', encoding='utf-8') as f:
    cluster_seq_id = pd.DataFrame(json.load(f)) 

try:
    os.makedirs(medaka_output)
except FileExistsError:
    shutil.rmtree(medaka_output)
    os.makedirs(medaka_output)

# Medaka paralelization
def process_file(file):
    seq_id = file.split(".")[0].split("/")[-1].split("_")[0]

    try:
        cluster_id = cluster_seq_id[cluster_seq_id["seq"] == seq_id].iloc[0]["cluster"]
    except IndexError:
        return f">{seq_id}|DUPCOUNT=0\n\n"

    cluster_id = cluster_id.split("_")[0]

    subprocess.run(
        [
            "medaka_consensus",
            "-g",
            "-d",
            f"{output_dir}/clusters_consensus/{cluster_id}.fasta",
            "-i",
            file,
            "-t",
            "1",
            "-o",
            f"{medaka_output}/{seq_id}",
        ]
    )

    try:
        dup_count = len(SeqIO.to_dict(SeqIO.parse(file, "fasta")))
    except ValueError as e:
        print(e)
        return None

    subprocess.run(
        [
            "medaka_variant",
            "-i",
            f"{output_dir}/clusters_consensus/{cluster_id}.fasta",
            "-r",
            f"{medaka_output}/{seq_id}/consensus.fasta",
            "-o",
            f"{medaka_output}/{seq_id}/variant_call/",
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    records = list(
        SeqIO.parse(f"{medaka_output}/{seq_id}/consensus.fasta", "fasta")
    )

    try:
        medaka_cons_seq = records[0].seq
    except IndexError as e:
        return f">{seq_id}|DUPCOUNT={dup_count}\n\n"
    
    return f">{seq_id}|DUPCOUNT={dup_count}\n{medaka_cons_seq}\n"


with concurrent.futures.ProcessPoolExecutor(max_workers=6) as executor:
    results = executor.map(
        process_file, [f.strip() for f in open(news_final_clusters,"r")]
    )

for result in results:
    if result:
        final_consensus.write(result)
        

final_consensus.close()
