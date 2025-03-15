import glob
import argparse
import pandas as pd

# Create argument parser
parser = argparse.ArgumentParser(description='Process some files.')

# Add arguments
parser.add_argument('-i','--input_dir',required=True, help='Input directory path')
parser.add_argument('-o','--output_dir',required=True, help='Output directory path')

# Parse arguments
args = parser.parse_args()

airr_files = glob.glob("%s/*clone-pass.tsv"%args.input_dir)

add_columns = {
    'sequence_id':[],
    'sample_id':[],
    'subject_id':[],
    'seq_len':[],
    'clone_size_count':[],
    'clone_size_freq':[]
}

for airr in airr_files:
    subject_id = airr.split("/")[-1].split("_atleast-2_clone-pass.tsv")[0]
    sample_id = subject_id
    airr = pd.read_csv(airr, sep="\t", low_memory=False)
    for clone in airr.iloc:
        add_columns["subject_id"].append(subject_id)
        add_columns["sample_id"].append(sample_id)
        add_columns["seq_len"].append(len(clone["sequence"]))

        clone_size_count = len(airr[airr["clone_id"]==clone["clone_id"]])
        add_columns["clone_size_count"].append(clone_size_count)
        add_columns["clone_size_freq"].append(clone_size_count/airr.shape[0])

        new_columns = clone["sequence_id"].split("|")[1:]
        for column in new_columns:
            column = column.split("=")
            try:
                add_columns[column[0]].append(column[1])
            except KeyError:
                add_columns[column[0]] = [column[1]]

        add_columns["sequence_id"].append("%s%s"%(clone["sequence_id"].split("|")[0],sample_id))

    try:
        merged_airr = pd.concat([merged_airr, airr], ignore_index=True)
    except:
        merged_airr = airr


for key, value in add_columns.items():
    print(f"Length of {key}: {len(value)}")

add_columns = pd.DataFrame(add_columns)
merged_airr = merged_airr.drop(columns="sequence_id")
final_airr = pd.concat([merged_airr, add_columns], axis=1)

print("Table size: %s x %s"%(final_airr.shape[0], final_airr.shape[1]))
final_airr = final_airr.rename(columns={'DUPCOUNT':'duplicate_count','CPRIMER':'cprimer'})
#final_airr = final_airr.rename(columns={'DUPCOUNT':'duplicate_count','CPRIMER':'cprimer'})
final_airr.to_csv("all_phage_airr.tsv", sep='\t', index=False)
