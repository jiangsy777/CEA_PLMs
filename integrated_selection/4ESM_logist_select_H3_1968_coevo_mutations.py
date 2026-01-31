# -*- coding: utf-8 -*-
import math
import pickle
import pandas as pd
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
# import Translator

#HA_of_H1_human_CA07_zhangsenlaoshi(different from NCBI strain)
# ref_seq_nt = 'ATGAAGGCAATACTAGTAGTTCTGCTATATACATTTGCAACCGCAAATGCAGACACATTATGTATAGGTTATCATGCGAACAATTCAACAGACACTGTAGACACAGTACTAGAAAAGAATGTAACAGTAACACACTCTGTTAACCTTCTAGAAGACAAGCATAACGGGAAACTATGCAAACTAAGAGGGGTAGCCCCATTGCATTTGGGTAAATGTAACATTGCTGGCTGGATCCTGGGAAATCCAGAGTGTGAATCACTCTCCACAGCAAGCTCATGGTCCTACATTGTGGAAACACCTAGTTCAGACAATGGAACGTGTTACCCAGGAGATTTCATCGATTATGAGGAGCTAAGAGAGCAATTGAGCTCAGTGTCATCATTTGAAAGGTTTGAGATATTCCCCAAGACAAGTTCATGGCCCAATCATGACTCGAACAAAGGTGTAACGGCAGCATGTCCTCATGCTGGAGCAAAAAGCTTCTACAAAAATTTAATATGGCTAGTTAAAAAAGGAAATTCATACCCAAAGCTCAGCAAATCCTACATTAATGATAAAGGGAAAGAAGTCCTCGTGCTATGGGGCATTCACCATCCACCTACTAGTGCTGACCAACAAAGTCTCTATCAGAATGCAGATGCATATGTTTTTGTGGGGTCATCAAGATACAGCAAGAAGTTCAAGCCGGAAATAGCAATAAGACCCAAAGTGAGGGATCAAGAAGGGAGAATGAACTATTACTGGACACTAGTAGAGCCGGGAGACAAAATAACATTCGAAGCAACTGGAAATCTAGTGGTACCGAGATATGCATTCGCAATGGAAAGAAATGCTGGATCTGGTATTATCATTTCAGATACACCAGTCCACGATTGCAATACAACTTGTCAAACACCCAAGGGTGCTATAAACACCAGCCTCCCATTTCAGAATATACATCCGATCACAATTGGAAAATGTCCAAAATATGTAAAAAGCACAAAATTGAGACTGGCCACAGGATTGAGGAATATCCCGTCTATTCAATCTAGAGGCCTATTTGGGGCCATTGCCGGTTTCATTGAAGGGGGGTGGACAGGGATGGTAGATGGATGGTACGGTTATCACCATCAAAATGAGCAGGGGTCAGGATATGCAGCCGACCTGAAGAGCACACAGAATGCCATTGACGAGATTACTAACAAAGTAAATTCTGTTATTGAAAAGATGAATACACAGTTCACAGCAGTAGGTAAAGAGTTCAACCACCTGGAAAAAAGAATAGAGAATTTAAATAAAAAAGTTGATGATGGTTTCCTGGACATTTGGACTTACAATGCCGAACTGTTGGTTCTATTGGAAAATGAAAGAACTTTGGACTACCACGATTCAAATGTGAAGAACTTATATGAAAAGGTAAGAAGCCAGCTAAAAAACAATGCCAAGGAAATTGGAAACGGCTGCTTTGAATTTTACCACAAATGCGATAACACGTGCATGGAAAGTGTCAAAAATGGGACTTATGACTACCCAAAATACTCAGAGGAAGCAAAATTAAACAGAGAAGAAATAGATGGGGTAAAGCTGGAATCAACAAGGATTTACCAGATTTTGGCGATCTATTCAACTGTCGCCAGTTCATTGGTACTGGTAGTCTCCCTGGGGGCAATCAGTTTCTGGATGTGCTCTAATGGGTCTCTACAGTGTAGAATATGTATTTAA'
# ref_seq = Translator.cds_translator(ref_seq_nt)
# print(ref_seq)

#HA_of_H3_human_PlosOne_H3_A_Aichi_2_1968_H3N2_aligned
ref_seq = 'MKTIIALSYIFCLALGQDLPGNDNSTATLCLGHHAVPNGTLVKTITDDQIEVTNATELVQSSSTGKICNN-PHRILDGIDCTLIDALLGDPHCDVFQNE-TWDLFVERSK-AFSNCYPYDVPDYASLRSLVASSGTLEFITEGF---TWTGVTQN-GGSNACKR-GPGSGFFSRLNWLTKSG--STYPVLNVTMPNNDNFDKLYIWGIHHPSTNQEQTSLYVQASGRVTVSTRRSQQTIIPNIGSRPWVRGLSSRISIYWTIVKPGDVLVINSNGNLIAPRGYFKMRT----------GKSSIMRSDAPIDTCISECITPNGSIPNDKPFQNVNKITYGACPKYVKQNTLKLATGMRNVPEKQT----RGLFGAIAGFIENGWEGMIDGWYGFRHQNSEGTGQAADLKSTQAAIDQINGKLNRVIEKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTRRQLRENAEDMGNGCFKIYHKCDNACIESIRNGTYDHDVYRDEALNNRFQIKGVELKSGYKDWILWISFAISCFLLCVVLLGFIMWACQRGNIRCNICI'
ref_seq = ref_seq.replace('-','~')
ref_seq_without_empty = ref_seq.replace('~','')
without_empty_site_lst = []
for i in range(len(ref_seq)):
    AA_site0 = ref_seq[i]
    if AA_site0 != '~':
        without_empty_site_lst.append(i + 1)
real_site_lst = [i+1 for i in range(len(ref_seq_without_empty))]
dict_ref_site_without_empty_site = dict(zip(without_empty_site_lst,real_site_lst))

# df_coevo = pd.read_csv('./generated_sequences/HA_H1_human_CA07_ref_site12_AA_pair_top5_site2_for_each_selected_site1.csv')
df_coevo = pd.read_csv('./generate/HA_H3_human_PlosOne_ref_site12_AA_pair_top3_site2_for_each_selected_site1.csv')

df_esm = pd.read_csv('../ESM-2_ft/PlosOne_H3_A_Aichi_2_1968_ref_HA_logist_new_virus10w_freeze32_lr3e-5/HA_mutation_effects_esm2_ft1800.csv')
esm_mutant_lst = list(df_esm['mutation'])
esm_score_lst = list(df_esm['delta_logit'])
dict_esm_mutant_score = dict(zip(esm_mutant_lst,esm_score_lst))

site1_ref_real_from1_lst = ['CA07_PlosOne_aligned_ref']
site2_ref_real_from1_lst = ['CA07_PlosOne_aligned_ref']
site1_mutant_score_lst = ['CA07_PlosOne_aligned_ref']
site2_mutant_score_lst = ['CA07_PlosOne_aligned_ref']
site12_mutant_score_lst = ['CA07_PlosOne_aligned_ref']


for  i in range(1,len(df_coevo)):
    site1_from1 = int(df_coevo.loc[i,'site1_from1_ref'])
    site2_from1 = int(df_coevo.loc[i, 'site2_from1_ref'])
    site1_from0_real = dict_ref_site_without_empty_site[site1_from1] - 1
    site2_from0_real = dict_ref_site_without_empty_site[site2_from1] - 1

    true_site1_aa_ref = ref_seq_without_empty[site1_from0_real]
    true_site2_aa_ref = ref_seq_without_empty[site2_from0_real]
    site1_mut = df_coevo.loc[i,'site1_mut']
    site1_mut_before = site1_mut[0]
    site1_mut_after = site1_mut[-1]
    if site1_mut_before != true_site1_aa_ref:
        print(site1_mut_before,true_site1_aa_ref)
        print(str(site1_from0_real),'wrong!!!')

    site1_mutant = site1_mut_before+str(site1_from0_real+1)+site1_mut_after
    site1_mutant_esm_score = float(dict_esm_mutant_score[site1_mutant])

    site2_mut = df_coevo.loc[i,'site2_mut']
    site2_mut_before = site2_mut[0]
    site2_mut_after = site2_mut[-1]
    if site2_mut_before != true_site2_aa_ref:
        print(site2_mut_before,true_site2_aa_ref)
        print(str(site2_from0_real),'wrong!!!')

    site2_mutant = site2_mut_before+str(site2_from0_real+1)+site2_mut_after
    site2_mutant_esm_score = float(dict_esm_mutant_score[site2_mutant])

    site1_ref_real_from1_lst.append(site1_from0_real + 1)
    site2_ref_real_from1_lst.append(site2_from0_real + 1)
    site1_mutant_score_lst.append(site1_mutant_esm_score)
    site2_mutant_score_lst.append(site2_mutant_esm_score)
    site12_mutant_score_lst.append(site1_mutant_esm_score+site2_mutant_esm_score)

df_coevo['site1_ref_real_from1'] = site1_ref_real_from1_lst
df_coevo['site2_ref_real_from1'] = site2_ref_real_from1_lst

df_coevo['site1_mutant_score'] = site1_mutant_score_lst
df_coevo['site2_mutant_score'] = site2_mutant_score_lst
df_coevo['site12_mutant_score'] = site12_mutant_score_lst
# df_coevo.to_csv('./generated_sequences/HA_H1_human_CA07_ref_site12_AA_pair_top5_site2_for_each_selected_site1_with_esm2_logist_score.csv',index=False)
df_coevo.to_csv('./generate/HA_H3_human_PlosOne_ref_site12_AA_pair_top3_site2_for_each_selected_site1_with_esm2_logist_score.csv',index=False)


