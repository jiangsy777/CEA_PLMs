# -*- coding: utf-8 -*-
import math
import pickle
import pandas as pd
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
# import Translator


#HA_of_H7_PlosOne(A_Netherlands_219_2003)
#PlosOne序列
# ref_seq = 'DKICLGHHAVSNGTKVNTLTERGVEVVNATETVERTNVPRICSK~GKRTVDLGQCGLLGTITGPPQCDQFLEF~SADLIIERRE~GSDVCYPGKFVNEEALRQILRESGGIDKETMGF~~~TYSGIRTN~GATSACRR~~SGSSFYAEMKWLLSNTDNAAFPQMTKSYKNTRKDPALIIWGIHHSGSTTEQTKLYGSGNKLITVGSSNYQQSFVPSPGARPQVNGQSGRIDFHWLMLNPNDTVTFSFNGAFIAPDRASFLR~~~~~~~~~~GKSMGIQSSVQVDANCEGDCYHSGGTIISNLPFQNINSRAVGKCPRYVKQESLMLATGMKNVPEIPKG~~~R'
#NCBI标准株
# ref_seq = 'DKICLGHHAVSNGTKVNTLTERGVEVVNATETVERTNVPRICSK~GKRTVDLGQCGLLGTITGPPQCDQFLEF~SADLIIERRE~GSDVCYPGKFVNEEALRQILRESGGIDKETMGF~~~TYSGIRTN~GTTSACRR~~SGSSFYAEMKWLLSNTDNAAFPQMTKSYKNTRKDPALIIWGIHHSGSTTEQTKLYGSGNKLITVGSSNYQQSFVPSPGARPQVNGQSGRIDFHWLILNPNDTVTFSFNGAFIAPDRASFLR~~~~~~~~~~GKSMGIQSEVQVDANCEGDCYHSGGTIISNLPFQNINSRAVGKCPRYVKQESLLLATGMKNVPEIPKR~~~R'
ref_seq = 'MNTQILVFALVASIPTNADKICLGHHAVSNGTKVNTLTERGVEVVNATETVERTNVPRICSK~GKRTVDLGQCGLLGTITGPPQCDQFLEF~SADLIIERRE~GSDVCYPGKFVNEEALRQILRESGGIDKETMGF~~~TYSGIRTN~GTTSACRR~~SGSSFYAEMKWLLSNTDNAAFPQMTKSYKNTRKDPALIIWGIHHSGSTTEQTKLYGSGNKLITVGSSNYQQSFVPSPGARPQVNGQSGRIDFHWLILNPNDTVTFSFNGAFIAPDRASFLR~~~~~~~~~~GKSMGIQSEVQVDANCEGDCYHSGGTIISNLPFQNINSRAVGKCPRYVKQESLLLATGMKNVPEIPKR~~~RRRGLFGAIAGFIENGWEGLIDGWYGFRHQNAQGEGTAADYKSTQSAIDQITGKLNRLIEKTNQQFELIDNEFTEVERQIGNVINWTRDSMTEVWSYNAELLVAMENQHTIDLADSEMNKLYERVKRQLRENAEEDGTGCFEIFHKCDDDCMASIRNNTYDHSKYREEAIQNRIQIDPVKLSSGYKDVILWFSFGASCFILLAIAMGLVFICVKNGNMRCTICI'
ref_seq = ref_seq.replace('-','~')
print(ref_seq.replace('~',''))
# a = ref_seq.replace('~','')
# for i in range(len(a)):
#     a1 = a[i]
#     a2 = ref0[i]
#     if a1!=a2:
#         print(i)
#         print(a1,a2)

################################
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

# 输入的fasta最后要回车一下

site1_lst = [118,123,126,116,122,132,148,147,117,149,193,189,190,228,217,224,225,226,222,213]
RBD_type_lst = ['130loop' for i in range(10)]+['190helix' for i in range(3)]+['220loop' for i in range(7)]
site1_RBD_sort_set = [i for i in range(1,11)]+[i for i in range(1,4)]+[i for i in range(1,8)]
start_loc = 19
# start_loc=1
#top_coevo_sites_num = 5
top_coevo_sites_num = 3
AA_pair_freq_thresh = 0.1
# site2_lst = [76,154,122,76,122,122,109,114,116,146,221,185,207,222,154,237,221,96,118,222,107]

file = '../data/generate/H13579_PlosOne_aligned_and_sampled_sequences.fasta'

fasta_file = parse_fasta(file)
seq_id_lst, seq_seq_lst = fasta_file[0], fasta_file[1]
# print(seq_id_lst)
# print(seq_seq_lst)
ref = str(seq_seq_lst[0])
# print('ref=',ref)
seqnum, colnum0 = len(seq_id_lst), len(seq_seq_lst[1])

colseq_lst = []
top_AA_pair_num_lst = []

for j in range(colnum0):
    colseq = []
    for k in range(seqnum):
        colseq.append(seq_seq_lst[k][j])
    colseq_lst.append(colseq)

site1_from1_plosone_lst = ['CA07_PlosOne_aligned_ref']
site2_from1_plosone_lst = ['CA07_PlosOne_aligned_ref']
site1_from1_ref_lst = ['CA07_PlosOne_aligned_ref']
site2_from1_ref_lst = ['CA07_PlosOne_aligned_ref']
site1_RBD_type_lst = ['CA07_PlosOne_aligned_ref']
site1_RBD_sort_lst = ['CA07_PlosOne_aligned_ref']
site2_sort_lst = ['CA07_PlosOne_aligned_ref']
site1_mut_lst = ['CA07_PlosOne_aligned_ref']
site2_mut_lst = ['CA07_PlosOne_aligned_ref']
mut_type_lst = ['CA07_PlosOne_aligned_ref']
mut_ref_seq_lst = [ref_seq]
mut_ref_seq_without_empty_lst = [ref_seq.replace('~','')]


for j in range(len(site1_lst)):
    site1 = site1_lst[j]-1
    RBD_type = RBD_type_lst[j]
    site1_RBD_sort = site1_RBD_sort_set[j]

    df_site1 = pd.read_csv('./generate/site_coevo_with_site'+str(site1+1)+'.csv')
    site2_lst = list(df_site1['site2_from1'])
    site2_lst.remove(site1+1)
    site1_ref_from0 = site1+start_loc-1
    colseq1 = colseq_lst[site1]
    for i in range(0,top_coevo_sites_num):
        site2_sort = i+1
        site2 = site2_lst[i]-1
        site2_ref_from0 = site2+start_loc-1
        colseq2 = colseq_lst[site2]
        dict_site12_AApair_count_all = {}
        for k in range(len(colseq1)):
            site12_AA_pair = str(colseq1[k])+str(colseq2[k])
            try:
                dict_site12_AApair_count_all[site12_AA_pair] += 1
            except:
                dict_site12_AApair_count_all[site12_AA_pair] = 1
        dict_site12_AApair_count_all = dict(sorted(dict_site12_AApair_count_all.items(),key=lambda x: x[1],reverse=True))
        dict_site12_AApair_count_top = {}
        for AApair0 in dict_site12_AApair_count_all.keys():
            count0 = dict_site12_AApair_count_all[AApair0]
            if count0 >= AA_pair_freq_thresh*seqnum:
                dict_site12_AApair_count_top[AApair0] = count0
        top_AA_pair_num_lst.append(len(dict_site12_AApair_count_top))

        ref_site1_AA = ref_seq[site1_ref_from0]
        ref_site2_AA = ref_seq[site2_ref_from0]
        ref_site12_AApair = str(ref_site1_AA)+str(ref_site2_AA)
        if ref_site12_AApair not in dict_site12_AApair_count_top.keys():
            print(site1+1,site2+1)
            print(site2_ref_from0+1)

        str_dict = str(dict_site12_AApair_count_top)
        if '~' not in str_dict:
            for AApair1 in dict_site12_AApair_count_top.keys():
                if (AApair1[0] != ref_site12_AApair[0])&(AApair1[1] != ref_site12_AApair[1]):#保证突变替换为双位点突变
                    AA_site1_pair = AApair1[0]
                    AA_site2_pair = AApair1[1]
                    replace_only_site1_site1_mut = ref_site1_AA+'_'+AA_site1_pair
                    replace_only_site1_site2_mut = ref_site2_AA + '_' + ref_site2_AA
                    ref_seq_only_site1_replace = ref_seq[:site1_ref_from0]+AA_site1_pair+ref_seq[site1_ref_from0+1:]
                    site1_from1_plosone_lst.append(site1+1)
                    site2_from1_plosone_lst.append(site2+1)
                    site1_from1_ref_lst.append(site1_ref_from0+1)
                    site2_from1_ref_lst.append(site2_ref_from0+1)
                    site1_RBD_type_lst.append(RBD_type)
                    site1_RBD_sort_lst.append(site1_RBD_sort)
                    site2_sort_lst.append(site2_sort)
                    site1_mut_lst.append(replace_only_site1_site1_mut)
                    site2_mut_lst.append(replace_only_site1_site2_mut)
                    mut_type_lst.append('site1_mut_only')
                    mut_ref_seq_lst.append(ref_seq_only_site1_replace)
                    if len(ref_seq_only_site1_replace.replace('~','')) != len(ref_seq.replace('~','')):
                        print('site1_mut_only_wrong')
                    mut_ref_seq_without_empty_lst.append(ref_seq_only_site1_replace.replace('~',''))

                    replace_only_site2_site1_mut = ref_site1_AA + '_' + ref_site1_AA
                    replace_only_site2_site2_mut = ref_site2_AA + '_' + AA_site2_pair
                    ref_seq_only_site2_replace = ref_seq[:site2_ref_from0] + AA_site2_pair + ref_seq[site2_ref_from0 + 1:]
                    site1_from1_plosone_lst.append(site1 + 1)
                    site2_from1_plosone_lst.append(site2 + 1)
                    site1_from1_ref_lst.append(site1_ref_from0 + 1)
                    site2_from1_ref_lst.append(site2_ref_from0+1)
                    site1_RBD_type_lst.append(RBD_type)
                    site1_RBD_sort_lst.append(site1_RBD_sort)
                    site2_sort_lst.append(site2_sort)
                    site1_mut_lst.append(replace_only_site2_site1_mut)
                    site2_mut_lst.append(replace_only_site2_site2_mut)
                    mut_type_lst.append('site2_mut_only')
                    mut_ref_seq_lst.append(ref_seq_only_site2_replace)
                    if len(ref_seq_only_site2_replace.replace('~','')) != len(ref_seq.replace('~','')):
                        print('site2_mut_only_wrong')
                    mut_ref_seq_without_empty_lst.append(ref_seq_only_site2_replace.replace('~', ''))

                    replace_site12_site1_mut = ref_site1_AA + '_' + AA_site1_pair
                    replace_site12_site2_mut = ref_site2_AA + '_' + AA_site2_pair
                    ref_seq_site12_replace = ref_seq_only_site1_replace[:site2_ref_from0] + AA_site2_pair + ref_seq_only_site1_replace[site2_ref_from0 + 1:]
                    site1_from1_plosone_lst.append(site1 + 1)
                    site2_from1_plosone_lst.append(site2 + 1)
                    site1_from1_ref_lst.append(site1_ref_from0 + 1)
                    site2_from1_ref_lst.append(site2_ref_from0+1)
                    site1_RBD_type_lst.append(RBD_type)
                    site1_RBD_sort_lst.append(site1_RBD_sort)
                    site2_sort_lst.append(site2_sort)
                    site1_mut_lst.append(replace_site12_site1_mut)
                    site2_mut_lst.append(replace_site12_site2_mut)
                    mut_type_lst.append('site12_mut_both')

                    mut_ref_seq_lst.append(ref_seq_site12_replace)
                    if len(ref_seq_site12_replace.replace('~','')) != len(ref_seq.replace('~','')):
                        print('site12_mut_both_wrong')
                    mut_ref_seq_without_empty_lst.append(ref_seq_site12_replace.replace('~', ''))
df_final = pd.DataFrame()
df_final['site1_from1_plosone'] = site1_from1_plosone_lst
df_final['site2_from1_plosone'] = site2_from1_plosone_lst
df_final['site1_from1_ref'] = site1_from1_ref_lst
df_final['site2_from1_ref'] = site2_from1_ref_lst
df_final['site1_RBD_type'] = site1_RBD_type_lst
df_final['site1_RBD_sort'] = site1_RBD_sort_lst
df_final['site2_sort'] = site2_sort_lst
df_final['site1_mut'] = site1_mut_lst
df_final['site2_mut'] = site2_mut_lst
df_final['mut_type'] = mut_type_lst
df_final['mut_ref_seq'] = mut_ref_seq_lst
df_final['mut_ref_seq_without_empty'] = mut_ref_seq_without_empty_lst
df_final.to_csv('./generate/HA_H7_avian_PlosOne_ref_site12_AA_pair_top'+str(top_coevo_sites_num)+'_site2_for_each_selected_site1.csv',index=False)

f1 = open('./generate/HA_H7_avian_PlosOne_ref_site12_AA_pair_top'+str(top_coevo_sites_num)+'_site2_for_each_selected_site1.fasta','w')

for i in range(len(df_final)):
    site1_from1_ref = df_final.loc[i,'site1_from1_ref']
    site2_from1_ref = df_final.loc[i,'site2_from1_ref']
    site1_mut = df_final.loc[i,'site1_mut']
    site2_mut = df_final.loc[i,'site2_mut']
    seq = df_final.loc[i,'mut_ref_seq']
    seq_id = str(site1_from1_ref)+'|'+str(site2_from1_ref)+'|'+site1_mut+'|'+site2_mut
    print('>'+seq_id,file=f1)
    print(seq,file=f1)
f1.close()
