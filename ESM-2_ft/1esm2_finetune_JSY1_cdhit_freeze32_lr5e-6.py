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

# 1. 模型和分词器加载
model_name = "./esm2_t33_650M_UR50D"
tokenizer = AutoTokenizer.from_pretrained(model_name)
model = AutoModelForMaskedLM.from_pretrained(model_name)
for layer in model.esm.encoder.layer[:32]:
    for param in layer.parameters():
        param.requires_grad = False
# 2. 数据准备

# 标准20种氨基酸字母表
amino_acids = "ACDEFGHIKLMNPQRSTVWY"

path = './'

sequences = []
'''
for file in os.listdir(path):
    if file.endswith('.csv') & file.startswith('seq'):
        genename = 'gene'
        print (genename)
        df = pd.read_csv(path + file)
        df = shuffle(df, random_state = 42)
        sequences = sequences + df[genename].tolist()
'''
for file in os.listdir(path):
    if file.endswith('.csv') & file.startswith('df_HA_with_aa_seq_H13579_human_avian_without_wrong_seq_cdhit999'):
        genename = 'HA_AA'
        print (genename)
        df = pd.read_csv(path + file)
        df = shuffle(df, random_state = 42)
        sequences = sequences + df[genename].tolist()

num_sequences = len(sequences)
print (num_sequences)

dataset = Dataset.from_dict({"sequence": sequences})

# # 3. 数据预处理函数

def preprocess_function(examples):
    # 分词（添加特殊 tokens，如<cls>和<<eos>）
    return tokenizer(
        examples["sequence"],
        padding="max_length",
        truncation=False,
        max_length=600,  # 根据模型支持的最大长度调整
        return_tensors="pt"
    )


tokenized_dataset = dataset.map(
    preprocess_function,
    batched=True,
    remove_columns=["sequence"]  # 移除原始序列列
)

# 4. 数据collator（用于MLM任务的掩码处理）
data_collator = DataCollatorForLanguageModeling(
    tokenizer=tokenizer,
    mlm=True,  # 启用掩码语言模型任务
    mlm_probability=0.15  # 随机掩码15%的token
)

# 5. 训练参数设置
training_args = TrainingArguments(
    output_dir="./esm2_finetuned_HA_H13579_CDHIT999_freeze32layers_lr5e-6",  # 模型保存路径
    overwrite_output_dir=True,
    num_train_epochs=10,  # 模拟数据量小， epochs可设小些
    per_device_train_batch_size=10,
    gradient_accumulation_steps=20,
    learning_rate=5e-6,
    logging_dir="./logs",
    logging_steps=50,  # 模拟数据少，日志步长设小些
    save_steps=200,
    fp16=True,  # 启用混合精度训练（需GPU支持）
    report_to="none"
)

# 6. 初始化Trainer并微调
trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=tokenized_dataset,
    data_collator=data_collator
)

# 开始微调
trainer.train()
