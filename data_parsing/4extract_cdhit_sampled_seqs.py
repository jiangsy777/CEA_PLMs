def parse_fasta(fasta_file):
    fr = open(fasta_file, 'r')
    contents = fr.readlines()
    idlst, seqlst, seq, seqnum = [], [], '', 0

    for line in contents:
        if line[:1] != '>':
            seq = seq + line[:-1]
        else:
            seqnum += 1
            seqlst.append(seq)
            seq = ''
            seqid = line[1:-1]
            idlst.append(seqid)
    seqlst = seqlst[1:]

    seqlast, line_num_all = '', len(contents)

    for j in range(line_num_all - 1, 0, -1):
        line = contents[j]
        if line[:1] != '>':
            seqlast = line[:-1] + seqlast
        else:
            break

    seqlst.append(seqlast)
    return idlst, seqlst

import pandas as pd

df_all = pd.read_csv('df_H13579_aligned_with_PlosOne_before_sampling.csv')
id_lst_all = list(df_all['seq_id'])
seq_lst_all = list(df_all['seq_HA'])
dict_id_seq = dict(zip(id_lst_all,seq_lst_all))
id_lst,seq_lst = [],[]

f1 = '../data/cdhit/H1_aligned_with_PlosOne_after_CDHIT_sampling_995.fasta'
file1 = parse_fasta(f1)
id_lst1,seq_lst1 = file1[0],file1[1]
for id0 in id_lst1:
    seq0 = dict_id_seq[id0]
    id_lst.append(id0)
    seq_lst.append(seq0.replace('-','~'))

f2 = '../data/cdhit/H3_aligned_with_PlosOne_after_CDHIT_sampling_995.fasta'
file2 = parse_fasta(f2)
id_lst2,seq_lst2 = file2[0],file2[1]
for id0 in id_lst2:
    seq0 = dict_id_seq[id0]
    id_lst.append(id0)
    seq_lst.append(seq0.replace('-','~'))

f3 = '../data/cdhit/H5_aligned_with_PlosOne_after_CDHIT_sampling_995.fasta'
file3 = parse_fasta(f3)
id_lst3,seq_lst3 = file3[0],file3[1]
for id0 in id_lst3:
    seq0 = dict_id_seq[id0]
    id_lst.append(id0)
    seq_lst.append(seq0.replace('-','~'))

f4 = '../data/generate/H135_cdhit_sampling/H7_aligned_with_PlosOne.fasta'
file4 = parse_fasta(f4)
id_lst4,seq_lst4 = file4[0],file4[1]
for id0 in id_lst4:
    seq0 = dict_id_seq[id0]
    id_lst.append(id0)
    seq_lst.append(seq0.replace('-','~'))

f5 = '../data/generate/H135_cdhit_sampling/H9_aligned_with_PlosOne.fasta'
file5 = parse_fasta(f5)
id_lst5,seq_lst5 = file5[0],file5[1]
for id0 in id_lst5:
    seq0 = dict_id_seq[id0]
    id_lst.append(id0)
    seq_lst.append(seq0.replace('-','~'))

f_all = open('../data/generate/H13579_PlosOne_aligned_and_sampled_sequences.fasta','w')
print(len(id_lst1),len(id_lst2),len(id_lst3),len(id_lst4),len(id_lst5))
print(len(id_lst),len(seq_lst))
for i in range(len(id_lst)):
    id1 = id_lst[i]
    seq1 = seq_lst[i]
    print('>'+id1,file=f_all)
    print(seq1,file=f_all)
f_all.close()
