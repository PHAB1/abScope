import pandas as pd
import glob
import json
import os
import argparse

parser = argparse.ArgumentParser(description="Pre Medaka script")
parser.add_argument('-j', '--json_path')
parser.add_argument('-l', '--list')
parser.add_argument('-c', '--chain')
parser.add_argument('-nfc', '--new_clusters')
parser.add_argument('-i', '--igcicle_out')
parser.add_argument('-jo', '--json_out')
args = parser.parse_args()

min_seqs_cluster = 3

# import snakemake inputs
jsonPath=args.json_path # get json file with id mapping
news_clusters_cdr3_done=args.list # done file, uc, get extention
igcicle_out=args.igcicle_out # igcicle output
json_out=args.json_out # igcicle output
igblast=pd.read_csv(igcicle_out)

# import snakemake params
chain_type=args.chain

# import snakemake outputs
news_final_clusters_dir=args.new_clusters
output_dir = "/".join(news_final_clusters_dir.split("/")[0:-3])

#cluster_seq_id = pd.DataFrame(cluster_seq_id)
with open(jsonPath, 'r', encoding='utf-8') as f:
    cluster_seq_id = pd.DataFrame(json.load(f))

dict_new_cluster_id = {
    'seq': [],
    'cdr_cluster_id': []
}

# get uc file from news cluster .fasta list (done.txt)
new_clusters_ucFiles = []
with open(news_clusters_cdr3_done, "r") as done: 
    for row in done:
        cleaned_row = row.strip("\n")
        base_name = os.path.splitext(cleaned_row)[0]
        new_clusters_ucFiles.append(".".join([base_name, "uc"]))

# Create final clusters dir
try:
    os.makedirs(news_final_clusters_dir)
except FileExistsError:
    pass

# Dictionary to get original draft for medaka consensus
#for file in glob.glob("%s/*.uc"%news_clusters_cdr3_dir):
for file in new_clusters_ucFiles:
    try:
        cluster = pd.read_csv(file,sep="\t",header=None)
    except:
        continue
    if cluster.empty:
        continue
    match_id_cluster = cluster[cluster.iloc[:, 0] == "H"].iloc[:,8:]
    filtered_cluster = cluster[cluster.iloc[:, 0] == "C"]
    filtered_cluster = filtered_cluster[filtered_cluster[2] > min_seqs_cluster]
    filtered_cluster = filtered_cluster.sort_values(by=2, ascending=False)
    print("Clusters Ordened By clusters lenght: \n%s"%filtered_cluster)
    for clus_id in filtered_cluster[8]:
        seqs_clus = match_id_cluster[match_id_cluster[9] == clus_id][8]
        dict_new_cluster_id["seq"].extend(seqs_clus)
        dict_new_cluster_id["cdr_cluster_id"].extend([clus_id]*len(seqs_clus))
        
        merge_clus_ig = pd.merge(seqs_clus, igblast, left_on=8, right_on='sequence_id', how='inner')
        output_cluster_file = open("%s/%s.fasta"%(news_final_clusters_dir,clus_id.split(";")[0]),"w")

        rec_id_dic = {}
        for rec in merge_clus_ig.iloc:

            if rec["sequence_id"] in rec_id_dic.keys():
                rec.id = "%s_%s"%(rec["sequence_id"], rec_id_dic[rec["sequence_id"]])
            else:
                rec.id = rec["sequence_id"]
                
            try:
                rec_id_dic[rec["sequence_id"]] += 1
            except:
                rec_id_dic[rec["sequence_id"]] = 1

            output_cluster_file.write(">%s\n%s\n"%(rec.id, rec["sequence"]))
        
        output_cluster_file.close()

df_new_cluster_id = pd.DataFrame(dict_new_cluster_id)
df_new_cluster_id = pd.merge(cluster_seq_id,df_new_cluster_id, left_on='seq', right_on='seq', how='left')
df_new_cluster_id.to_csv("%s/meta_%s.csv"%(output_dir,chain_type), index=False)

with open(json_out, 'w', encoding='utf-8') as f:
    # Converter o DataFrame para uma lista de dicionários (orient='records')
    data_to_dump = df_new_cluster_id.to_dict(orient='records')

    # Usar json.dump() para escrever no arquivo com ensure_ascii=False e indent
    json.dump(data_to_dump, f, ensure_ascii=False, indent=4)
