import Translator
import pandas as pd

df_all = pd.read_csv('../data/df_IAV_HA_deduplicated_labels_98787.csv')
serotype_H_lst = []
for i in range(len(df_all)):
    print(i)
    serotype0 = str(df_all.loc[i,'Serotype'])
    host0 = df_all.loc[i,'Host1']
    if 'H' in serotype0:
        serotype_lst0 = serotype0.split('N')
        serotype_H = serotype_lst0[0]
    else:
        serotype_H = 'other'
    serotype_H_lst.append(serotype_H)
df_all['Serotype_H'] = serotype_H_lst
df_H1 = df_all[(df_all['Serotype_H']=='H1')&(df_all['Host1']=='Human')]
df_H3 = df_all[(df_all['Serotype_H']=='H3')&(df_all['Host1']=='Human')]
df_H5 = df_all[(df_all['Serotype_H']=='H5')&(df_all['Host1']=='Avian')]
df_H7 = df_all[(df_all['Serotype_H']=='H7')&(df_all['Host1']=='Avian')]
df_H9 = df_all[(df_all['Serotype_H']=='H9')&(df_all['Host1']=='Avian')]
df_sele = pd.concat([df_H1,df_H3,df_H5,df_H7,df_H9],ignore_index=True)
aa_seq_lst = []
flag_lst = []
for i in range(len(df_sele)):
    print(i)
    nt_seq = df_sele.loc[i,'HA']
    aa_seq0 = Translator.cds_translator(nt_seq)
    if aa_seq0[-1]=='*':
        aa_seq = aa_seq0[:-1]
    else:
        aa_seq = aa_seq0
    aa_seq_lst.append(aa_seq)
    if ('~' not in aa_seq)&('*' not in aa_seq)&(aa_seq!='It is not a cds sequence!'):
        flag_lst.append('T')
    else:
        flag_lst.append('F')
df_sele['HA_AA'] = aa_seq_lst
df_sele['flag_seq_quality'] = flag_lst

df_sele_final = df_sele[df_sele['flag_seq_quality']=='T']
df_sele_final.to_csv('../data/generate/df_HA_with_aa_seq_H13579_human_avian_without_wrong_seq.csv',index=False)

a = 'MENIVLLLAIVSLVKSDQICIGYHANNSTEQVDTIMEKNVTVTHAQDILEKTHNGKLCDLNGVKPLILKDCSVAGWLLGNPMCDEFIRVPEWSYIVERANPANDLCYPGSLNDYEELKHMLSRINHFEKIQIIPKSSWPNHETSLGVSAACPYQGAPSFFRNVVWLIKKNDAYPTIKISYNNTNREDLLILWGIHHSNNAEEQTNLYKNPITYISVGTSTLNQRLAPKIATRSQVNGQRGRMDFFWTILKPDDAIHFESNGNFIAPEYAYKIVKKGDSTIMKSGVEYGHCNTKCQTPVGAINSSMPFHNIHPLTIGECPKYVKSNKLVLATGLRNSPLREKRRKRGLFGAXAGFIEGGWQGMVDGWYGYHHSNEQGSGYAADKESTQKAIDGVTNKVNSIIDKMNTQFEAVGREFNNLERRIENLNKKMEDGFLDVWTYNAELLVLMENERTLDFHDSNVKNLYDKVRLQLRDNAKELGNGCFEFYHKCDNECMESVRNGTYDYPQYSEEARLKREEISGVKLESVGTYQILSIYSTAASSLALAIMMAGLSLWMCSNGSLQCRICI'

df_sele_final0 = pd.read_csv('../data/generate/df_HA_with_aa_seq_H13579_human_avian_without_wrong_seq.csv')
sero_set1 = ['H1','H3']
for seroh in sero_set1:
    f1 = open('../data/generate/HA_with_aa_seq_without_wrong_seq_'+seroh+'.fasta','w')
    df_sele_final = df_sele_final0[df_sele_final0['Serotype_H']==seroh]
    df_sele_final.to_csv('../data/generate/df_HA_with_aa_seq_'+seroh+'_without_wrong_seq.csv',index=False)
    df_sele_final = pd.read_csv('../data/generate/df_HA_with_aa_seq_'+seroh+'_without_wrong_seq.csv')
    print('>A_Texas_37_2024_H5N1_HA',file=f1)
    print(a,file=f1)
    for i in range(len(df_sele_final)):
        id = df_sele_final.loc[i,'strain_name']
        sero_H = df_sele_final.loc[i,'Serotype_H']
        host = df_sele_final.loc[i,'Host1']
        seq = df_sele_final.loc[i,'HA_AA']
        seq_id = id+'_'+host+'_'+sero_H
        print('>'+seq_id,file=f1)
        print(seq,file=f1)
    f1.close()

f1 = open('../data/generate/HA_with_aa_seq_without_wrong_seq_H579.fasta','w')
df_sele_final = df_sele_final0[(df_sele_final0['Serotype_H']=='H5')|(df_sele_final0['Serotype_H']=='H7')|(df_sele_final0['Serotype_H']=='H9')]
df_sele_final.to_csv('../data/generate/df_HA_with_aa_seq_H579_without_wrong_seq.csv',index=False)
df_sele_final = pd.read_csv('../data/generate/df_HA_with_aa_seq_H579_without_wrong_seq.csv')
print('>A_Texas_37_2024_H5N1_HA',file=f1)
print(a,file=f1)
for i in range(len(df_sele_final)):
    id = df_sele_final.loc[i,'strain_name']
    sero_H = df_sele_final.loc[i,'Serotype_H']
    host = df_sele_final.loc[i,'Host1']
    seq = df_sele_final.loc[i,'HA_AA']
    seq_id = id+'_'+host+'_'+sero_H
    print('>'+seq_id,file=f1)
    print(seq,file=f1)
f1.close()

