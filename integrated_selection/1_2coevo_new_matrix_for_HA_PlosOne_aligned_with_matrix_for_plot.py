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


def high_occur_pair_visualize(path,sort_num,site1,site2,seqnum,RBS_name):
    import matplotlib.pyplot as plt
    import os

    colseq1 = colseq_lst[site1]
    colseq2 = colseq_lst[site2]
    thresh_freq = 0.01
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

    dict_Ccodonpairs = dict(sorted(dict_Ccodonpairs.items(),key=lambda x: x[1],reverse=True))
    print('dict_Ccodonpairs=', dict_Ccodonpairs)
    pair_number0 = len(dict_Ccodonpairs)
    if pair_number0 == 1:
        coevo_score = 1
        for pair0, count0 in dict_Ccodonpairs.items():
            pair_ref = pair0
    else:
        max_count = max(dict_Ccodonpairs.values())
        for pair0, count0 in dict_Ccodonpairs.items():
            if count0 == max_count:
                pair_ref = pair0
                continue
    title = 'pair_ref: '+pair_ref+'_'+str(site1+1)+'_'+str(site2+1)+'_RBS'+str(RBS_name)
    title1 = 'sort_num'+str(sort_num)+'_'+'site2_from1'+str(site2 + 1) + '_RBS' + str(RBS_name)
    dict_high_occur_pair_freq = {}
    for pair000, count000 in dict_Ccodonpairs.items():
        freq000 = count000 / seqnum
        if freq000 > thresh_freq:
            dict_high_occur_pair_freq[pair000] = count000
    names = list(dict_high_occur_pair_freq.keys())
    values = list(dict_high_occur_pair_freq.values())
    plt.figure(figsize=(6,4))
    plt.bar(names,values)
    plt.title(title,fontsize=14)
    plt.tight_layout()
    folder_path = path+'result/fig/RBS'+str(RBS_name)
    if not os.path.exists(folder_path):
        os.makedirs(folder_path,exist_ok=True)
    file_format = "png"
    plt.savefig(f"{folder_path}/{title1}.{file_format}",dpi=300)

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
    print('dict_Ccodonpairs1111111111=', dict_Ccodonpairs)
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

# RBS_name_lst = [i for i in range(115,133)]
# RBS_name_lst = [i for i in range(133,150)]
# RBS_name_lst = [i for i in range(188,196)]
# RBS_name_lst = [i for i in range(213,230)]
# RBS_name_lst = [237,248]
# RBS_name_lst = [190]
start_loc = 70  # 70 (10AA before 80)-1(python from 0)
end_loc = 290  # 290 (10AA after 280)
RBS_name_lst = [i for i in range(start_loc-1,end_loc)]

site_10_lst = [i for i in RBS_name_lst]
coevo_matrix = []
for lll in range(len(RBS_name_lst)):
    RBS_name = RBS_name_lst[lll]
    site_10 = site_10_lst[lll]
    site_1 = site_10
    site_1_and_other_coevo_score_lst = []
    start_loc = 70 #70 (10AA before 80)-1(python from 0)
    end_loc = 290#290 (10AA after 280)
    for site_2 in range(start_loc-1,end_loc):
        sites12 = (str(site_1), str(site_2))
        colseq1 = colseq_lst[site_1]
        colseq2 = colseq_lst[site_2]
        codon_result0 = coevo_matrix_calculation(colseq1, colseq2, seqnum)
        codon_result = codon_result0[0]
        codon_score_result = codon_result[1]
        site_1_and_other_coevo_score_lst.append(codon_score_result)
    coevo_matrix.append(site_1_and_other_coevo_score_lst)

coevo_matrix = np.array(coevo_matrix)
np.save('./generate/coevo_matrix_H13579_PlosOne_70_290.npy',coevo_matrix)

    #
    # df_r = pd.DataFrame()
    # site_1_and_other_coevo_score_lst0 = [np.log10(c+1e-100) for c in site_1_and_other_coevo_score_lst]
    # if max(site_1_and_other_coevo_score_lst0)-min(site_1_and_other_coevo_score_lst0)>0:
    #     site_1_and_other_coevo_score_lst1 = [(c-min(site_1_and_other_coevo_score_lst0))/(max(site_1_and_other_coevo_score_lst0)-min(site_1_and_other_coevo_score_lst0)) for c in site_1_and_other_coevo_score_lst0]
    # else:
    #     site_1_and_other_coevo_score_lst1 = site_1_and_other_coevo_score_lst0
    # df_r['site1_from1'] = [RBS_name for a in range(len(site_1_and_other_coevo_score_lst))]
    # df_r['site2_from1'] = [a for a in range(start_loc,end_loc+1)]
    # df_r['site_coevo'] = site_1_and_other_coevo_score_lst1
    # df_r1 = df_r.sort_values(by='site_coevo', ascending=False)
    # df_r1.to_csv('./'+path+'result/site_coevo_with_site'+str(RBS_name)+'.csv', index=False)
    #

    #
    #
    # site2_lst = list(df_r1['site2_from1'])
    #
    # for i in range(10):
    #     site2_0 = site2_lst[i]
    #     print(site2_0)
    #     sort_num = i
    #     high_occur_pair_visualize(path,sort_num,site_1, site2_0-1, seqnum, RBS_name)
    #
    # str_len = 5
    # site2_from1_lst = list(df_r['site2_from1'])
    # coevo_lst = list(df_r['site_coevo'])
    # coevo_mean_lst = []
    # slide_start_site_from1_lst = []
    # for i in range(len(df_r)-str_len+1):
    #     slide_window_lst = coevo_lst[i:i+str_len]
    #     coevo_mean = np.mean(slide_window_lst)
    #     slide_start_site_from1 = site2_from1_lst[i]
    #     coevo_mean_lst.append(coevo_mean)
    #     slide_start_site_from1_lst.append(slide_start_site_from1)
    # df_r2 = pd.DataFrame()
    # df_r2['coevo_mean'] = coevo_mean_lst
    # df_r2['site_from1_slide_start'] = slide_start_site_from1_lst
    # df_r2 = df_r2.sort_values(by='coevo_mean',ascending=False)
    # df_r2.to_csv('./'+path+'result/sliding_window_site_coevo_mean_with_site'+str(RBS_name)+'.csv', index=False)
    #
    # plt.figure(figsize=(8,5))
    # plt.plot(slide_start_site_from1_lst,coevo_mean_lst,marker='o',linewidth=2,color='#1f77b4')
    #
    # plt.tight_layout()
    # plt.savefig('./'+path+'result/fig/RBS'+str(RBS_name)+'/sliding_window'+str(str_len)+'.png',dpi=300)
