import pandas as pd

df = pd.read_csv('Predicted_prob_H13579_PlosOne_aligned_cdhit98_test_seq_split_random30.csv')
f1 = open('Predicted_prob_H13579_PlosOne_aligned_cdhit98_test_seq_split_random30_sampled15_each.fasta','w')

sero_lst = []

for i in range(len(df)):
    seq_id1 = df.loc[i,'seq_id']
    seq_id_lst1 = seq_id1.split('_')
    sero_lst.append(seq_id_lst1[-1])
df['serotype'] = sero_lst

df1 = df[df['serotype'] == 'H1'].sample(n=15,random_state=42)
df3 = df[df['serotype'] == 'H3'].sample(n=15,random_state=42)
df5 = df[df['serotype'] == 'H5'].sample(n=15,random_state=42)
df7 = df[df['serotype'] == 'H7'].sample(n=15,random_state=42)
df9 = df[df['serotype'] == 'H9'].sample(n=15,random_state=42)

df_sampled = pd.concat([df1,df3,df5,df7,df9],ignore_index=True)


for i in range(len(df_sampled)):
    seq_id = df_sampled.loc[i,'seq_id']
    seq = df_sampled.loc[i,'seq']
    print('>'+seq_id,file=f1)
    print(seq.replace('~','-'),file=f1)
f1.close()