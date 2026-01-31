# -*- coding: utf-8 -*-
import math
import pickle
import pandas as pd
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt


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


################################
def coevo_matrix_calculation(colseq1, colseq2, seqnum):
    thresh_freq = 0.01
    thresh_freq_high = 0.05
    colseq_12 = []
    colseq_1 = []
    colseq_2 = []
    dict_Csite1 = {}
    dict_Csite2 = {}
    for i in range(len(colseq1)):
        pair_12 = colseq1[i] + colseq2[i]
        colseq_12.append(pair_12)
        colseq_1.append(colseq1[i])
        colseq_2.append(colseq2[i])
    dict_Ccodonpairs = {}
    for codonpair in colseq_12:
        try:
            dict_Ccodonpairs[codonpair] += 1
        except:
            dict_Ccodonpairs[codonpair] = 1
    print('dict_Ccodonpairs=', dict_Ccodonpairs)
    pair_number0 = len(dict_Ccodonpairs)
    if pair_number0 == 1:
        coevo_score = 1
    else:
        max_count = max(dict_Ccodonpairs.values())
        for pair0, count0 in dict_Ccodonpairs.items():
            if count0 == max_count:
                pair_ref = pair0
                continue
        two_sites_mutation_freq = 0
        for pair00 in dict_Ccodonpairs.keys():
            freq00 = dict_Ccodonpairs[pair00] / seqnum
            if (pair00[0] != pair_ref[0]) & (pair00[1] != pair_ref[1]) & (freq00 > thresh_freq):
                two_sites_mutation_freq += freq00

        if two_sites_mutation_freq == 0:
            coevo_score = 0
        else:
            site_pair_entropy = 0
            for pair000, count000 in dict_Ccodonpairs.items():
                freq000 = count000 / seqnum
                site_pair_entropy -= freq000 * math.log2(freq000)

            dict_Csite1 = {}
            for col1 in colseq_1:
                try:
                    dict_Csite1[col1] += 1
                except:
                    dict_Csite1[col1] = 1

            dict_Csite2 = {}
            for col2 in colseq_2:
                try:
                    dict_Csite2[col2] += 1
                except:
                    dict_Csite2[col2] = 1

            co_mut_score = 1
            high_occur_pair_num = 0
            high_occur_pair_num1 = 0
            for pair000, count000 in dict_Ccodonpairs.items():
                freq000 = count000 / seqnum
                high_occur_pair_num1 += 1
                if freq000 > thresh_freq:
                    high_occur_pair_num += 1
                if freq000 > thresh_freq_high:
                    count_col1 = dict_Csite1[pair000[0]]
                    count_col2 = dict_Csite2[pair000[1]]
                    co_mut_score = co_mut_score * ((count000 ** 2) / (count_col1 * count_col2))
                    print(co_mut_score,pair000,count000,count_col1,count_col2)
            # coevo_score = (two_sites_mutation_freq*co_mut_score)/site_pair_entropy
            coevo_score = (two_sites_mutation_freq**2 * co_mut_score**5) / (math.sqrt(2**site_pair_entropy)*high_occur_pair_num*math.log2(high_occur_pair_num1))
            print('two_sites_mutation_freq',two_sites_mutation_freq,'co_mut_score',co_mut_score,'site_pair_entropy',site_pair_entropy,'high_occur_pair_num',high_occur_pair_num,'high_occur_pair_num1',high_occur_pair_num1)
    print('coevo_score=', coevo_score)
    print('dict1=', dict_Csite1)
    print('dict2=', dict_Csite2)
    sitepair = (site_1, site_2)
    print('------------------------', site_1, '_', site_2, 'co_evo_score calculation end', '------------------------')
    site_score = (sitepair, coevo_score)
    # dict_Ccodonpairs_sorted = dict(sorted(dict_Ccodonpairs.items(), key=lambda e: e[1], reverse=True))
    # print('dict_Ccodonpairs=', dict_Ccodonpairs)
    # print('dict_Ccodonpairs_sorted=', dict_Ccodonpairs_sorted)
    # dict_Ccodonpairs_sorted_head = dict(list(dict_Ccodonpairs_sorted.items())[0:11])
    # print('dict_Ccodonpairs_sorted_head=', dict_Ccodonpairs_sorted_head)

    return site_score, dict_Ccodonpairs


# 输入的fasta最后要回车一下

file = '../data/generate/H13579_PlosOne_aligned_and_sampled_sequences.fasta'
path = ''

# file = './seq_H13579_70_280/sampled_H5_HA_RBD_for_coevo.fasta'
# path = 'H5_single_coevo/'

# codon
fasta_file = parse_fasta(file)
seq_id_lst, seq_seq_lst = fasta_file[0], fasta_file[1]
ref = str(seq_seq_lst[0])
# print('ref=',ref)
seqnum, colnum0 = len(seq_id_lst), len(seq_seq_lst[1])

colseq_lst = []
for j in range(colnum0):
    colseq = ''
    for k in range(seqnum):
        colseq += seq_seq_lst[k][j]
    colseq_lst.append(colseq)

seqlst = fasta_file[1]
seqlst_1st = seqlst[0]
seq_len = len(seqlst_1st)
sites = []
codon_co_evo_score = []
dict_site2_dictC = {}
dict_site12_dictC = {}
matrix_all_coevo_score = []

start_loc = 70  # 70(10AA before 80)-1(python from 0)
end_loc = 290  # 290(10AA after 280)
site1_lst = []
coevo_mean_lst = []
for site_1 in range(start_loc-1,end_loc):
    site_1_and_other_coevo_score_lst = []
    # site1_lst.append(site_1-12+1)
    # site1_lst.append(site_1)
    for site_2 in range(start_loc-1,end_loc):
        sites12 = (str(site_1), str(site_2))
        colseq1 = colseq_lst[site_1]
        colseq2 = colseq_lst[site_2]
        codon_result0 = coevo_matrix_calculation(colseq1, colseq2, seqnum)
        codon_result = codon_result0[0]
        codon_score_result = codon_result[1]
        site_1_and_other_coevo_score_lst.append(codon_score_result)
        # site_1_and_other_coevo_score_lst1 = [np.log10(c + 1e-100) for c in site_1_and_other_coevo_score_lst]
    # if site_1 == 260:
    #     print(site_1_and_other_coevo_score_lst)
    # colseq_set_num = len(list(set(colseq1)))
    if (1 not in site_1_and_other_coevo_score_lst)&(colseq1.count('~')<seqnum*0.05):
        site1_lst.append(site_1 + 1)
        top_lst_all = list(sorted(site_1_and_other_coevo_score_lst,reverse=True))
        top_lst = top_lst_all[1:]

        # site_1_and_other_coevo_score_lst1 = [np.log10(c+1e-100) for c in site_1_and_other_coevo_score_lst]
        coevo_mean_lst.append(np.mean(top_lst))
        #site_1_and_other_coevo_score_lst1 = [np.log10(c + 1e-100) for c in site_1_and_other_coevo_score_lst]
df = pd.DataFrame()

# df['site2'] = [kkk for kkk in range(70,281)]
# df['lst_196'] = a196
# df['lst_145'] = a145
# df.to_csv('145_196.csv',index=False)
df_r = pd.DataFrame()
# site_1_and_other_coevo_score_lst1 = [(c-min(site_1_and_other_coevo_score_lst0))/(max(site_1_and_other_coevo_score_lst0)-min(site_1_and_other_coevo_score_lst0)) for c in site_1_and_other_coevo_score_lst0]
df_r['site1_from1'] = site1_lst
df_r['coevo_mean_all_sites'] = coevo_mean_lst
df_r = df_r.sort_values(by='coevo_mean_all_sites',ascending=False)
df_r.to_csv('./'+path+'generate/selected_site_coevo.csv', index=False)
