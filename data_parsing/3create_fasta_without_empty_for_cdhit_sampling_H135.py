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

start_loc = 70 #80位点（头部第一个beta-sheet起始位约80（79）），向前约10AAs
end_loc = 290 #280位点（头部最后一个beta-sheet起始位约280），向后约10AAs，两个位点中间含空格

df_all = pd.DataFrame()
f1 = '../data/generate/HA_with_aa_seq_without_wrong_seq_H1_mafft_del_empty_with_plosone_final.fasta'
file1 = parse_fasta(f1)
id_lst1,seq_lst1 = file1[0],file1[1]
id_lst1,seq_lst1 = id_lst1[1:],seq_lst1[1:]

f2 = '../data/generate/HA_with_aa_seq_without_wrong_seq_H3_mafft_del_empty_with_plosone_final.fasta'
file2 = parse_fasta(f2)
id_lst2,seq_lst2 = file2[0],file2[1]
id_lst2,seq_lst2 = id_lst2[1:],seq_lst2[1:]

f3 = '../data/generate/HA_with_aa_seq_without_wrong_seq_H5_mafft_del_empty_with_plosone_final.fasta'
file3 = parse_fasta(f3)
id_lst3,seq_lst3 = file3[0],file3[1]
id_lst3,seq_lst3 = id_lst3[1:],seq_lst3[1:]

f4 = '../data/generate/HA_with_aa_seq_without_wrong_seq_H7_mafft_del_empty_with_plosone_final.fasta'
file4 = parse_fasta(f4)
id_lst4,seq_lst4 = file4[0],file4[1]
id_lst4,seq_lst4 = id_lst4[1:],seq_lst4[1:]

f5 = '../data/generate/HA_with_aa_seq_without_wrong_seq_H9_mafft_del_empty_with_plosone_final.fasta'
file5 = parse_fasta(f5)
id_lst5,seq_lst5 = file5[0],file5[1]
id_lst5,seq_lst5 = id_lst5[1:],seq_lst5[1:]
# print(len(id_lst5),len(set(id_lst5)))

id_lst = id_lst1+id_lst2+id_lst3+id_lst4+id_lst5
seq_lst = seq_lst1+seq_lst2+seq_lst3+seq_lst4+seq_lst5
seq_cut_lst = [seq[start_loc-1:end_loc] for seq in seq_lst]
seq_cut_without_empty_lst = [seq.replace('-','') for seq in seq_cut_lst]

df_all['seq_id'] = id_lst
df_all['seq_HA'] = seq_lst
df_all['seq_'+str(start_loc)+'_to_'+str(end_loc)] = seq_cut_lst
df_all['seq_'+str(start_loc)+'_to_'+str(end_loc)+'_without_empty'] = seq_cut_without_empty_lst
df_all.to_csv('df_H13579_aligned_with_PlosOne_before_sampling.csv',index=False)

fw1 = open('../data/generate/H135_cdhit_sampling/H1_aligned_with_PlosOne_before_sampling.fasta','w')
for i in range(len(id_lst1)):
    id0 = id_lst1[i]
    seq_complete0 = seq_lst1[i]
    seq_cut0 = seq_complete0[start_loc-1:end_loc]
    seq_cut_without_empty0 = seq_cut0.replace('-','')
    print('>'+id0,file=fw1)
    print(seq_cut_without_empty0,file=fw1)
fw1.close()

fw2 = open('../data/generate/H135_cdhit_sampling/H3_aligned_with_PlosOne_before_sampling.fasta','w')
for i in range(len(id_lst2)):
    id0 = id_lst2[i]
    seq_complete0 = seq_lst2[i]
    seq_cut0 = seq_complete0[start_loc-1:end_loc]
    seq_cut_without_empty0 = seq_cut0.replace('-','')
    print('>'+id0,file=fw2)
    print(seq_cut_without_empty0,file=fw2)
fw2.close()

fw3 = open('../data/generate/H135_cdhit_sampling/H5_aligned_with_PlosOne_before_sampling.fasta','w')
for i in range(len(id_lst3)):
    id0 = id_lst3[i]
    seq_complete0 = seq_lst3[i]
    seq_cut0 = seq_complete0[start_loc-1:end_loc]
    seq_cut_without_empty0 = seq_cut0.replace('-','')
    print('>'+id0,file=fw3)
    print(seq_cut_without_empty0,file=fw3)
fw3.close()

fw4 = open('../data/generate/H135_cdhit_sampling/H7_aligned_with_PlosOne.fasta','w')
for i in range(len(id_lst4)):
    id0 = id_lst4[i]
    seq_complete0 = seq_lst4[i]
    seq_cut0 = seq_complete0[start_loc-1:end_loc]
    seq_cut_without_empty0 = seq_cut0.replace('-','')
    print('>'+id0,file=fw4)
    print(seq_cut_without_empty0,file=fw4)
fw4.close()

fw5 = open('../data/generate/H135_cdhit_sampling/H9_aligned_with_PlosOne.fasta','w')
for i in range(len(id_lst5)):
    id0 = id_lst5[i]
    seq_complete0 = seq_lst5[i]
    seq_cut0 = seq_complete0[start_loc-1:end_loc]
    seq_cut_without_empty0 = seq_cut0.replace('-','')
    print('>'+id0,file=fw5)
    print(seq_cut_without_empty0,file=fw5)
fw5.close()