import pandas as pd
from Bio import SeqIO
import numpy as np
import random as rd
import subprocess
import shutil
import time
import sys
import re
import os

# Function for convert fasta files to data frame
def fasta2df(fasta_file):
    records = []
    for rec in SeqIO.parse(fasta_file, "fasta"):
        product=rec.description.split("|")[3]
        if product == "P" or product == "F":
            records.append({
                "id": (rec.name).split("|")[1],
                "sequence": str(rec.seq)
            })

    df = pd.DataFrame(records)
    return(df)

def df2fasta(df, out_file):
    with open(out_file, 'w') as f:
        for index, row in df.iterrows(): # Usando iterrows para acessar colunas por nome
            seq_id = row['sequence_id']
            sequence = row['sequence']
            f.write(f">{seq_id}\n")
            f.write(f"{sequence}\n")
    f.close()

# function to compile data frames from each V, D or J fasta file
def compileDf(fasta_files):
    df_list = []
    for fa in fasta_files:
        df_list.append(fasta2df(fa))
    return(pd.concat(df_list))

# function to create one aleatory VH V-D-J and VK V-J sequence
def recombine(df,ch_types):
    # make recombination, consider n sequences for each recombination #
    recombined_seq = {}
    for k,values in ch_types.items(): # for chain type
        combs = [k+v for v in values]
        rearr = []
        id_list = []
        for comb in combs:
            df_filtered = df[df['id'].str.startswith(comb)].copy() # possible rows for pattern (V,D or J)
            random_index = np.random.choice(df_filtered.index)
            frag = df_filtered.loc[random_index] # aleatory V,D or J fragments
            id_list.append(frag["id"]) # get V-D-J family
            rearr.append(frag["sequence"].replace(".","")) # get fragment in linear format (no align dots)
        id = "|".join(id_list)
        rearr = "".join(rearr)
        recombined_seq[id] = rearr
    return(recombined_seq)

# Run Igblastn - complete steps (Get dictionary and return AIRR format file)
def run_igblastn(dicScfv, ig_db_path, threads):
    # write fasta file for igblastn #
    input_fasta = "igblastn_input.fasta"

    # open file
    with open(input_fasta, "w") as fasta_file:
        # write sequences
        for scfv in dicScfv:
            [fasta_file.write(">%s\n%s\n" % (k, v)) for k, v in scfv.items()]
    # close file
    fasta_file.close()

    # run igblastn AIRR tabular annotated file using fasta file #
    subprocess.run(["igblastn", 
        "-germline_db_J", "%s/IGHLKJ_edit.fasta"%ig_db_path, 
        "-germline_db_V", "%s/IGHLKV_edit.fasta"%ig_db_path, 
        "-germline_db_D", "%s/IGHD_edit.fasta"%ig_db_path,
        "-domain_system", "kabat",
        "-query", "%s"%input_fasta, 
        "-outfmt", "19", # outfmt 19 for default table 
        "-show_translation", 
        "-auxiliary_data", "%s/human_gl.aux"%ig_db_path, 
        "-out", "ig_out.tsv",
        "-num_threads", "%s"%threads
        ])

    # read AIRR file
    airr_df = pd.read_csv("ig_out.tsv", sep="\t", header=0)
    return(airr_df)

# get the choosen nucleotide and return another one
# for example: if A is choosen, return T or C or G
def AID_like(nuc2mute,nucs):
    if nuc2mute.upper() in nucs:
        nucs.remove(nuc2mute)
        return rd.choice(nucs)
    else:
        return "Y"

# function to alter sequences in acordance with error or mutation
# fragmentos serão passados no caso de região especifico
def mutate(sequence,e=0,n_mut=0): 
    try:
        indices = list(range(len(sequence)))
    except TypeError:
        return("")
    if e > 0:
        randomIndex = rd.sample(indices, round(len(sequence) * e))
    elif n_mut > 0:
        randomIndex = rd.sample(indices, n_mut)
    for i in randomIndex:
        sequence = sequence[:(i-1)] + AID_like(nuc2mute=sequence[i-1],nucs=["A","T","C","G"]) + sequence[i:]
    return(sequence)

def mutBunch(ctrl_airr, lib_airr, replicates, clones, f_error_list, f_mut_list, n_mut_list, region="sequence"):
    #--> generalized <--#
    # error
    for e in f_error_list:
        for c in range(clones):
            error_rec=ctrl_airr.copy()
            error_rec["original_sequence"] = ctrl_airr["sequence"] # get original sequence
            recWerror_seq=[mutate(ctrl_airr.iloc[i][region],e=e) for i in range(len(ctrl_airr))]
            error_rec[region]=recWerror_seq
            error_rec["sequence_id"]=["%s|error=%s|clone=%s"%(id,e,c) for id in error_rec["sequence_id"]]
            lib_airr = pd.concat([lib_airr,error_rec], ignore_index=True)

    # Mutation
    for n_mut in n_mut_list: # para cada n mutantes presentes na seq
        mut_rec=ctrl_airr.copy()
        mut_rec["original_sequence"] = ctrl_airr["sequence"] # get original sequence
        ctrl_sequences = [ctrl_airr.iloc[i][region] for i in range(len(ctrl_airr))] # para cada airr row captura sequencia
        recWmut_seq=[mutate(seq,e=0,n_mut=n_mut) for seq in ctrl_sequences] # mutate each sequence 
        mut_rec[region]=recWmut_seq
        mut_rec["sequence_id"]=["%s|mut=%s"%(id,n_mut) for id in mut_rec["sequence_id"]]
        lib_airr = pd.concat([lib_airr,mut_rec], ignore_index=True)
    
    return(lib_airr)

# function to generate aleatory recombined V-D-J library
# considering: freq error, specific mutant freq (for n mutants)
# generalized or Fw and CDR cases in both illumina/nanopore paths
def gen_lib(df,ch_types,replicates,clones,f_error_list,f_mut_list,n_mut_list,ig_db_path):
    # 8 cases including generalized and region specific
    #ab_generated = 
    #print("Thre will be generated %s sequences"%ab_generated)

    regions_list = [
        "sequence",
        "fwr1",
        "cdr1",
        "fwr2",
        "cdr2",
        "fwr3",
        "cdr3",
        "fwr4"]

    # Create here the recombined sequences
    control_rec=[recombine(df, ch_types) for i in range(replicates)]
    
    #--> Create here the igblastn AIRR tabular annotated file and read <--#
    # get AIRR format annotated sequences #
    ctrl_airr = run_igblastn(control_rec, ig_db_path, 2)

    #--> for all sequence and Fw1/2/3/4 or CDR1/2/3 region specific - Error prone and mutation simulations <--#
    lib_airr_list = []
    for region in regions_list:
        lib_airr = mutBunch(ctrl_airr, ctrl_airr.copy(), replicates, clones, f_error_list, f_mut_list, n_mut_list, region=region) 
        if region !="sequence":
            lib_airr["sequence_id"] = lib_airr["sequence_id"].astype(str) + "|region=" + region
            lib_airr["sequence"] = lib_airr['fwr1'].astype(str) + lib_airr['cdr1'].astype(str) + lib_airr['fwr2'].astype(str) + lib_airr['cdr2'].astype(str) + lib_airr['fwr3'].astype(str) + lib_airr['cdr3'].astype(str) + lib_airr['fwr4'].astype(str)
        lib_airr_list.append(lib_airr)
    
    lib_airr = pd.concat(lib_airr_list, ignore_index=True)

    # For mutated rows create identical sequences clone times substituting ids eg.: |clones=n
    # Use -> for row with 'mut=' -> create new table with the same sequence and diff; id ('|clone') and merge with the original table
    # Lista para armazenar as novas linhas
    new_rows = []

    # Itera pelas linhas com "mut=" no sequence_id
    for _, row in lib_airr[lib_airr['sequence_id'].str.contains("mut=")].iterrows():
        for f_mut in f_mut_list:
            n_mut = round(f_mut * clones)
            n_orig = clones - n_mut  # garante que total == clones

            # Sequências mutadas
            for i in range(1, n_mut + 1):
                mut_row = row.copy()
                mut_row['sequence'] = row['sequence']
                mut_row['original_sequence'] = row['original_sequence']
                mut_row['sequence_id'] = f"{row['sequence_id']}|clone={i}|f_mut={f_mut}|original=FALSE"
                new_rows.append(mut_row)

            # Sequências originais
            for i in range(1, n_orig + 1):
                orig_row = row.copy()
                orig_row['sequence'] = row['original_sequence']
                orig_row['original_sequence'] = row['original_sequence']
                orig_row['sequence_id'] = f"{row['sequence_id']}|clone={i}|f_mut={f_mut}|original=TRUE"
                new_rows.append(orig_row)

    # Cria o novo DataFrame com as linhas expandidas
    lib_expanded_airr = pd.DataFrame(new_rows)
    lib_expanded_airr = pd.concat([lib_airr, lib_expanded_airr], ignore_index=True)

    # return final lib
    return(lib_expanded_airr)

# generate nanopore and illumina samples 
# "GGTGGTGGTGGTTCTGGTGGTGGTGGTCTGGTGGTGGTGGTTCT" as linker for scFv
# cases of error or mutants plus control (2 samples for each run)
# 2 file for nanopore or 4 files for illumina run (treatment and control)
def generate_sample(lib, pattern, out_dir, illumina=True):
    ## control and tratment ##
    # filter by pattern for VH
    treat_vh = lib[(lib['sequence_id'].str.contains(pattern)) & (lib['locus'].str.contains("IGH", na=False))].reset_index(drop=True)
    ctrl_vh = lib[(lib['original_sequence'].fillna('').str.strip() == '') & (lib['locus'].str.contains("IGH", na=False))].reset_index(drop=True)

    # filter by pattern for VL
    treat_vl = lib[(lib['sequence_id'].str.contains(pattern)) & (lib['locus'].str.contains("IGK", na=False))].reset_index(drop=True)
    ctrl_vl = lib[(lib['original_sequence'].fillna('').str.strip() == '') & (lib['locus'].str.contains("IGK", na=False))].reset_index(drop=True)

    # mkdir 
    try:
        os.makedirs(out_dir)
    except FileExistsError:
        shutil.rmtree(out_dir)
        os.makedirs(out_dir)
    
    pattern_name = pattern.split("=")[0]
    if illumina:
        # write fasta samples treatment and control VH and VH (4 files)
        # treatment
        df2fasta(treat_vh, out_file=f"{out_dir}/VH_illu_treat_{pattern_name}.fasta")
        df2fasta(treat_vl, out_file=f"{out_dir}/VL_illu_treat_{pattern_name}.fasta")
        # control
        df2fasta(ctrl_vh, out_file=f"{out_dir}/VH_illu_ctrl_{pattern_name}.fasta")
        df2fasta(ctrl_vl, out_file=f"{out_dir}/VL_illu_ctrl_{pattern_name}.fasta")

    # nanopore
    else:
        # set linker
        linker = "GGTGGTGGTGGTTCTGGTGGTGGTGGTCTGGTGGTGGTGGTTCT"
        
        # join using "§" for id and linker for sequence #
        # treatment
        treat_merged = treat_vh.copy()
        treat_merged["sequence_id"] = treat_merged["sequence_id"] + '§' + treat_vl['sequence_id']
        treat_merged['sequence'] = treat_merged['sequence'] + linker + treat_vl['sequence']

        # control
        ctrl_merged = ctrl_vh.copy()
        ctrl_merged["sequence_id"] = ctrl_merged["sequence_id"] + '§' + ctrl_vl['sequence_id']
        ctrl_merged['sequence'] = ctrl_merged['sequence'] + linker + ctrl_vl['sequence']

        # write fasta samples treatment 
        df2fasta(treat_merged, out_file=f"{out_dir}/nano_treat_{pattern_name}.fasta")

        # write fasta samples control
        df2fasta(ctrl_merged, out_file=f"{out_dir}/nano_ctrl_{pattern_name}.fasta")

def main():
    # check if ig_db_path is provided
    if len(sys.argv) < 2:
        print("Uso: python3 simulateVDJ.py <ig_db_path>")
        sys.exit(1)

    # get ig_db_path from command line argument
    ig_db_path = sys.argv[1]

    # Lib generation parmeters
    replicates = 10 # number of V-D-J recombinations
    clones = 5 # number of clones for replicate (each replicate is V-D-J recombination, each clone is related to an error or mutation condition in their respective list)
    f_error_list = [0.01,0.05,0.1,0.15,0.2]
    f_mut_list = [0.01,0.05,0.1,0.15,0.2]
    n_mut_list = [1,2,3,4,5]

    # path to fasta files with V, D or J or V and J for light chain
    fasta_files = [f"{ig_db_path}/IGHV.fasta",
              f"{ig_db_path}/IGHD.fasta",
              f"{ig_db_path}/IGHJ.fasta",
              f"{ig_db_path}/IGKV.fasta",
              f"{ig_db_path}/IGKJ.fasta"]

    df = compileDf(fasta_files)

    ###--> aleatory recombination <--###
    # patterns associated in order  #
    # (dont change the order here ) #
    ch_types = {
           "IGH": ["V",
                   "D",
                   "J"],
           "IGK": ["V",
                   "J"]
           }

    # generate the library
    lib_expanded_airr = gen_lib(df,ch_types,replicates,clones,f_error_list,f_mut_list,n_mut_list,ig_db_path)

    # write lib_airr in tsv format
    lib_expanded_airr.to_csv("lib_airr.tsv", sep="\t", index=False)

    #--> export illumina (VH|VH separated and nanopore (VH|VL linked in scFv)) <--#
    # Create control sample with no "mut=" and no "error=" to comparation, for both seq. generations
    # for each seq. generation create an error prone sample and mutated sample (by selecting 'mut=' or 'error=')
    # Illumina #
    generate_sample(lib_expanded_airr, pattern="mut=", out_dir="samples/illumina/mut")
    generate_sample(lib_expanded_airr, pattern="error=", out_dir="samples/illumina/error")

    # nanopore #
    generate_sample(lib_expanded_airr, pattern="mut=", out_dir="samples/nano/mut", illumina=False)
    generate_sample(lib_expanded_airr, pattern="error=", out_dir="samples/nano/error", illumina=False)


if __name__ == "__main__":
    main()

