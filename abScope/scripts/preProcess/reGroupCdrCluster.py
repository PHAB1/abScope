import pandas as pd
import glob
import time

cluster_seq_id = pd.DataFrame(cluster_seq_id)

time.sleep(1)
dict_new_cluster_id = {
    'seq': [],
    'cdr_cluster_id': []
}

# Dictionary to get original draft for medaka consensus
for file in glob.glob("%s/*.uc"%news_clusters_cdr3_dir):
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
df_new_cluster_id.to_csv("%s/meta_%s.csv"%(output_dir,args.chain_type), index=False)