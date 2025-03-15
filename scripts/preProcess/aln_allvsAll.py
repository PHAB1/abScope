import pandas as pd
from Bio.Seq import Seq
from Bio import pairwise2

teste = pd.read_csv("quality-passbarcode10_VL_quality-pass.paf",sep="\t",header=None)
print(teste)
aln_dic = {}
for line in teste.iloc:
    id1 = line[0]
    id2 = line[5]
    match = int(line[9])
    aln_len = int(line[10])
    #t_len = int(line[6])
    if aln_len > 300:
        #print("%s, %s, %s, %s\n"%(line, match/(aln_len), match, aln_len))
        try:
            aln_dic[id1].append(id2)
        except KeyError:
            try:
                aln_dic[id2].append(id1)
            except KeyError:
                aln_dic[id1] = [id2]

list_count_len = []
for i in aln_dic.items():
    list_count_len.append(int(i))
    if int(len(i[1])) > 10:
        print("%s - count: %s"%(i[0],len(i[1])))
print(max(list_count_len))
