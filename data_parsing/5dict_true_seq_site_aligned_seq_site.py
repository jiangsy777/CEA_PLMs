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

f1 = '../data/generate/H13579_PlosOne_aligned_and_sampled_sequences.fasta'
file1 = parse_fasta(f1)
id_lst, seq_lst = file1[0],file1[1]

for i in range(len(id_lst)):
    id0 = id_lst[i]
    seq0 = seq_lst[i]

    if id0 == 'A_California_04_2009_H1N1':
        seq1 = seq0
        n1 = 17
        true_site_lst1, aligned_seq_lst1 = [], []
        for j in range(len(seq1)):
            if seq1[j] != '~':
                n1 += 1
                true_site_lst1.append(n1)
                aligned_seq_lst1.append(j + 1)
        df1 = pd.DataFrame()
        df1['true_site'] = true_site_lst1
        df1['aligned_site'] = aligned_seq_lst1
        df1.to_csv('../data/generate/H1_dict_true_aligned_site.csv', index=False)

    elif id0 == 'A_Aichi_2_1968_H3N2':
        seq1 = seq0
        n1 = 26
        true_site_lst1, aligned_seq_lst1 = [], []
        for j in range(len(seq1)):
            if seq1[j] != '~':
                n1 += 1
                true_site_lst1.append(n1)
                aligned_seq_lst1.append(j + 1)
        df1 = pd.DataFrame()
        df1['true_site'] = true_site_lst1
        df1['aligned_site'] = aligned_seq_lst1
        df1.to_csv('../data/generate/H3_dict_true_aligned_site.csv', index=False)

    elif id0 == 'A_Texas_37_2024_H5N1_HA':
        seq1 = seq0
        n1 = 4#alphafold预测序列从13开始
        true_site_lst1, aligned_seq_lst1 = [], []
        for j in range(len(seq1)):
            if seq1[j] != '~':
                n1 += 1
                true_site_lst1.append(n1)
                aligned_seq_lst1.append(j + 1)
        df1 = pd.DataFrame()
        df1['true_site'] = true_site_lst1
        df1['aligned_site'] = aligned_seq_lst1
        df1.to_csv('../data/generate/H5_dict_true_aligned_site.csv', index=False)

    elif id0 == 'A_Netherlands_219_2003':
        # seq1 = seq0
        #NCBI标准株
        seq1 = 'DKICLGHHAVSNGTKVNTLTERGVEVVNATETVERTNVPRICSK~GKRTVDLGQCGLLGTITGPPQCDQFLEF~SADLIIERRE~GSDVCYPGKFVNEEALRQILRESGGIDKETMGF~~~TYSGIRTN~GTTSACRR~~SGSSFYAEMKWLLSNTDNAAFPQMTKSYKNTRKDPALIIWGIHHSGSTTEQTKLYGSGNKLITVGSSNYQQSFVPSPGARPQVNGQSGRIDFHWLILNPNDTVTFSFNGAFIAPDRASFLR~~~~~~~~~~GKSMGIQSEVQVDANCEGDCYHSGGTIISNLPFQNINSRAVGKCPRYVKQESLLLATGMKNVPEIPKR~~~R'
        n1 = 0
        true_site_lst1, aligned_seq_lst1 = [], []
        for j in range(len(seq1)):
            if seq1[j] != '~':
                n1 += 1
                true_site_lst1.append(n1)
                aligned_seq_lst1.append(j + 1)
        df1 = pd.DataFrame()
        df1['true_site'] = true_site_lst1
        df1['aligned_site'] = aligned_seq_lst1
        df1.to_csv('../data/generate/H7_dict_true_aligned_site.csv', index=False)

    elif id0 == 'A_Swine_HK_9_1998_H9N2':
        seq1 = seq0
        n1 = 0
        true_site_lst1, aligned_seq_lst1 = [], []
        for j in range(len(seq1)):
            if seq1[j] != '~':
                n1 += 1
                true_site_lst1.append(n1)
                aligned_seq_lst1.append(j + 1)
        df1 = pd.DataFrame()
        df1['true_site'] = true_site_lst1
        df1['aligned_site'] = aligned_seq_lst1
        df1.to_csv('../data/generate/H9_dict_true_aligned_site.csv', index=False)






