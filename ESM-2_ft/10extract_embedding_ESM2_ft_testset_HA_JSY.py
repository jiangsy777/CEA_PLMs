import pandas as pd
import os
from sklearn.utils import shuffle
import string
import random
import torch
from datasets import Dataset
from transformers import (
    AutoModelForMaskedLM,
    AutoTokenizer,
    TrainingArguments,
    Trainer,
    DataCollatorForLanguageModeling
)
import torch
import Translator
import pickle

 # 1. 从微调后的模型提取序列Embedding
def get_sequence_embedding(sequence, model, tokenizer):
    """输入单个蛋白序列，返回其Embedding（使用<cls> token的输出）"""
    model.eval()
    with torch.no_grad():
        inputs = tokenizer(
            sequence,
            #padding="max_length",
            truncation=True,
            #max_length=max_length,
            return_tensors="pt"
        ).to(device)  # 移动到模型所在设备（GPU/CPU）

        outputs = model(**inputs, output_hidden_states=True)
        # 取最后一层的<cls> token输出作为序列Embedding
        cls_embedding = outputs.hidden_states[-1][:, :, :].squeeze()
        cls_embedding = cls_embedding[1:-1]
    return cls_embedding.cpu().numpy()  # 转为numpy数组


# extract_embedding
model_name = "./esm2_finetuned_virus10w_freeze32layers_lr3e-5/checkpoint-1800" 
#model_name = "./esm2_t33_650M_UR50D"
tokenizer = AutoTokenizer.from_pretrained(model_name)
model = AutoModelForMaskedLM.from_pretrained(model_name)
cds_lst = ['HA']
device = "cuda" if torch.cuda.is_available() else "cpu"
model = model.to(device)
for cds in cds_lst:
    pkl_path = './'+cds+'_sampled_Serotype_new_not_in_training_set_and_1wHA_ems2_ft_embedding.pkl'
    embed_lst = []
    print(cds+'_finished!')
    df_test = pd.read_csv ('./'+cds+'_sampled_Serotype_new_not_in_training_set_and_1wHA.csv')
    for i in range(len(df_test)):
        print(i)
        test_sequence1 = str(df_test.loc[i,'CDS_4_amino_seq'])
        #test_sequence1 = Translator.cds_trans(test_sequence0)
        if '~' in test_sequence1:
            test_sequence1 = test_sequence1.replace('~','')
        #test_sequence = test_sequence1[:-1]
        test_sequence = test_sequence1
        print(len(test_sequence))
        embedding = get_sequence_embedding(
            test_sequence,
            model=model,
            tokenizer=tokenizer
        )
        embed_lst.append(embedding)
    with open(pkl_path,'wb') as f_pkl:
        pickle.dump(embed_lst,f_pkl)   
    print(cds+'_finished!!!!!!')     
#print(f"测试序列: {embedding[-1]}...（长度100）")
#print(f"序列Embedding形状: {embedding.shape}")  # 输出应为(384,)（对应8M模型）
