import subprocess
import pandas as pd
from Bio.SeqIO import FastaIO
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import shutil
import time
import os
import argparse

from Bio import SeqIO

parser = argparse.ArgumentParser(description="Process fasta file.")
parser.add_argument('file', type=str, help="The path to the file to be processed")
args = parser.parse_args()

ig_db_path = "references/database/human/ig_db"
input_fasta = args.file
threads = 12
output_dir = "results/igCicle/%s"%input_fasta.split("/")[-1].split(".")[0]

positionsAb_end = ["cdr2_end", 
                    "fwr3_end", 
                    "cdr3_end",
                    "fwr4_end"]

#positionsAb_end = ["fwr4_end"] 

try:
    os.makedirs(output_dir)
except FileExistsError:
    shutil.rmtree(output_dir)
    os.makedirs(output_dir)

# Separate VH and VK sequences
vh_file = open("%s/VH_seqs.fasta"%output_dir,"w")
vl_file = open("%s/VL_seqs.fasta"%output_dir,"w")

start_time = time.time()
count_ab = 1
while count_ab < 3:
    print("Starting Cicle %s"%count_ab)
    # Run igBlastn
    print("Runing igblast")
    subprocess.run(["igblastn", 
        "-germline_db_J", "%s/IGHLKJ_edit.fasta"%ig_db_path, 
        "-germline_db_V", "%s/IGHLKV_edit.fasta"%ig_db_path, 
        "-germline_db_D", "%s/IGHD_edit.fasta"%ig_db_path,
        "-domain_system", "kabat",
        "-query", "%s"%input_fasta, 
        "-outfmt", "19", # outfmt 19 for default table 
        "-show_translation", 
        "-auxiliary_data", "%s/human_gl.aux"%ig_db_path, 
        "-out", "%s/%s_igCicle.tsv"%(output_dir,count_ab),
        "-num_threads", "%s"%threads
        ])
    print("End of igblast: %s sec"%(time.time()-start_time))

    cicleAnnot = pd.read_csv("%s/%s_igCicle.tsv"%(output_dir,count_ab), delimiter='\t')

    try:
        meta_Annot_table = pd.merge(meta_Annot_table, cicleAnnot, on='sequence_id', how='outer')
    except:
        meta_Annot_table = cicleAnnot

    seq_to_new_cicle = open("%s/cicle_%s.fasta"%(output_dir, count_ab), "w")
    passed_seq_cicle = open("%s/cicle_%s_passed.fasta"%(output_dir, count_ab), "w")
    
    print("Starting Cicle Annotation")
    for rec in SeqIO.parse(input_fasta, "fasta"):
        try:
            Annot_seq = (cicleAnnot[cicleAnnot["sequence_id"] == "%s"%(rec.id)]) # add id plus cicle in header
        except KeyError:
            print("Not found ID:", rec.id)
        
        Annot_seq_end_pos = Annot_seq.filter(items=positionsAb_end).dropna(axis=1)
        
        Annot_seq_start_pos = Annot_seq["fwr1_start"]
        Annot_seq_end_pos = Annot_seq_end_pos.max().max()

        rec.seq = Seq("".join(str(Annot_seq.iloc[0]["sequence"]).split("-"))) # "germline_alignment"? , but "sequence" normaly

        if (not np.isnan(Annot_seq_end_pos)) and not Annot_seq_start_pos.isna().any():
            start = int(Annot_seq_start_pos.iloc[0])
            end = int(Annot_seq_end_pos)
            seq_copy = Seq(str(rec.seq))
            
            # Write Rec To Atual cicle - Fasta Passed sequences
            passed_seq_cicle.write(">%s\n%s\n"%(rec.id,rec.seq[int(start)-0:int(end)+180]))
            
            # Write Rec To Next cicle - Fasta
            rec.seq = seq_copy
            seq_to_new_cicle.write(">%s\n%s\n"%(rec.id,rec.seq[:int(start)] + rec.seq[int(end):]))

            # Separate VH from VK
            if Annot_seq.iloc[0]["locus"] == "IGH":
                vh_file.write(">%s\n%s\n"%(rec.id, rec.seq[int(start)-0:int(end)+180]))
            else:
                vl_file.write(">%s\n%s\n"%(rec.id, rec.seq[int(start)-0:int(end)+180]))

            # annotated meta only for passed sequences
            Annot_seq = Annot_seq.iloc[0]

    print("Eding Cicle annotation")
    passed_seq_cicle.close()
    seq_to_new_cicle.close()

    input_fasta = "%s/cicle_%s.fasta"%(output_dir,count_ab)
    
    print("End of Cicle %s: %s sec"%(count_ab, (time.time()-start_time)))
    count_ab += 1

# Close VH and VK/VL files
vh_file.close()
vl_file.close()

# Def To Filter meta columns based on the pattern
def colInPattern(col):
    pattern = ["sequence_id",
            "locus",
            "v_call",
            "j_call",
            "cdr3"]

    for p in pattern:
        if p in col:
            return True
    return False

print("Starting meta annotation")
filtered_meta_cols = []
filtered_meta_cols.append([str(col) for col in meta_Annot_table.columns if colInPattern(col)])

# Create a new DataFrame with only the filtered pattern columns names
meta_Annot_table = meta_Annot_table[filtered_meta_cols[0]]

# conserve only first allele in j and v call
# cicle 1
meta_Annot_table['j_call_x'] = meta_Annot_table['j_call_x'].str.split('*').str[0]
meta_Annot_table['v_call_x'] = meta_Annot_table['v_call_x'].str.split('*').str[0]
# cicle 2
meta_Annot_table['j_call_y'] = meta_Annot_table['j_call_y'].str.split('*').str[0]
meta_Annot_table['v_call_y'] = meta_Annot_table['v_call_y'].str.split('*').str[0]

## Create meta file ##
meta_Annot_table["ig_c1"] = (meta_Annot_table["locus_x"] + "_" +
                            meta_Annot_table["j_call_x"] + "_" +
                            meta_Annot_table["v_call_x"]
                            )

meta_Annot_table["ig_c2"] = (meta_Annot_table["locus_y"] + "_" +
                            meta_Annot_table["j_call_y"] + "_" +
                            meta_Annot_table["v_call_y"]
                            )

meta_Annot_table = meta_Annot_table[["sequence_id", "ig_c1", "ig_c2"]]
meta_Annot_table.to_csv("%s/meta.tsv"%output_dir,sep="\t",index=False)

print("Total time: %s"%(time.time()-start_time))
