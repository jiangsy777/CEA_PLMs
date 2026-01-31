import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import confusion_matrix, classification_report, roc_curve, auc
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
from torch.nn.functional import softmax
import warnings
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore")

# 设置全局字体为Arial
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False


# ===================== 1. 核心工具函数 =====================
def parse_label_from_seq_id(seq_id):
    """从seq_id解析标签: xxx_human_H5N5 → 1; xxx_avian_H5N1 → 0"""
    parts = seq_id.split("_")
    label_str = parts[-2].lower()
    label_str1 = parts[-1]
    return 1 if (label_str == "human") | (label_str1 in ['H1N1', 'H3N2']) else 0


def pad_embeddings(embedding_list, pad_value=0.0):
    """补齐变长embedding（词数×768）到统一长度，展平为一维向量"""
    max_word_num = max([emb.shape[0] for emb in embedding_list])  # 最大分词数
    padded_embeddings = []
    for emb in embedding_list:
        pad_len = max_word_num - emb.shape[0]
        # 仅在分词维度补齐，768维特征保持不变
        padded_emb = np.pad(emb, ((0, pad_len), (0, 0)), mode="constant", constant_values=pad_value)
        padded_embeddings.append(padded_emb.flatten())  # 展平：(词数×768,)
    return np.array(padded_embeddings)


# ===================== 2. 自定义数据集类 =====================
class SeqDataset(Dataset):
    """极简数据集：加载补齐+展平后的embedding和标签"""

    def __init__(self, embeddings, labels):
        self.embeddings = torch.tensor(embeddings, dtype=torch.float32)
        self.labels = torch.tensor(labels, dtype=torch.long)

    def __len__(self):
        return len(self.embeddings)

    def __getitem__(self, idx):
        return self.embeddings[idx], self.labels[idx]


# ===================== 3. 超轻量线性分类器 =====================
class UltraSimpleLinearClassifier(nn.Module):
    """仅 Dropout + 单层Linear，适配展平后的embedding维度"""

    def __init__(self, embedding_dim, dropout_rate=0.1):
        super().__init__()
        self.dropout = nn.Dropout(dropout_rate)
        self.fc = nn.Linear(embedding_dim, 2)  # 二分类：avian=0, human=1
        self.fc_features = None  # 存储全连接层输出，用于风险量化

    def forward(self, x):
        x = self.dropout(x)
        self.fc_features = self.fc(x)
        return self.fc_features


# ===================== 4. 带loss记录的快速训练函数 =====================
def fast_train_with_loss_record(model, train_loader, criterion, optimizer, epochs=20, device="cpu"):
    """训练模型并记录每个epoch的平均loss"""
    model.train()
    epoch_loss_list = []  # 存储每个epoch的平均loss
    for epoch in range(epochs):
        total_loss = 0.0
        for embeds, labels in train_loader:
            embeds, labels = embeds.to(device), labels.to(device)
            optimizer.zero_grad()
            outputs = model(embeds)
            loss = criterion(outputs, labels)
            loss.backward()
            optimizer.step()
            total_loss += loss.item()
        # 计算当前epoch的平均loss
        avg_loss = total_loss / len(train_loader)
        epoch_loss_list.append(avg_loss)
        print(f"Epoch [{epoch + 1}/{epochs}], Average Loss: {avg_loss:.4f}")
    # 保存loss到CSV：两列分别为epoch序号(从0开始)和对应loss
    loss_df = pd.DataFrame({
        "epoch": list(range(epochs)),
        "avg_loss": epoch_loss_list
    })
    loss_df.to_csv("epoch_loss_record.csv", index=False)
    print("Epoch loss记录已保存至 epoch_loss_record.csv")
    return model, epoch_loss_list


# ===================== 5. 数据加载与预处理 =====================
def load_and_preprocess_data(csv_path, npy_path, seq_id_col, test_size=0.3, random_state=42):
    """
    整合数据加载逻辑：
    1. CSV读取seq_id；2. NPY读取embedding；3. 一一匹配；4. 标签解析+embedding补齐
    """
    # 1. 读取CSV文件（仅含seq_id，无embedding）
    df = pd.read_csv(csv_path)
    # 2. 读取npy文件中的embedding列表（顺序需与CSV中seq_id完全一致！）
    embeddings = np.load(npy_path, allow_pickle=True)
    # 校验样本数量是否一致
    assert len(df) == len(embeddings), f"CSV样本数({len(df)})与NPY embedding数({len(embeddings)})不匹配！"

    # 3. 解析标签
    df["label"] = df[seq_id_col].apply(parse_label_from_seq_id)
    # 4. 补齐embedding并展平
    padded_embeds = pad_embeddings(embeddings)

    # 5. 分层划分训练集/测试集（保证类别分布一致）
    train_df, test_df = train_test_split(
        df, test_size=test_size, random_state=random_state, stratify=df["label"]
    )
    # 6. 提取对应训练集和测试集的embedding
    X_train = padded_embeds[train_df.index]
    y_train = train_df["label"].values
    X_test = padded_embeds[test_df.index]
    y_test = test_df["label"].values
    test_seq_ids = test_df[seq_id_col].values
    test_seqs = test_df["seq"].values
    return X_train, y_train, X_test, y_test, test_seq_ids, test_seqs


# ===================== 6. 训练集下采样平衡 =====================
def downsample_balance(X, y, random_state=42):
    """快速下采样：直接截断多数类，平衡训练集"""
    X_0, y_0 = X[y == 0], y[y == 0]
    X_1, y_1 = X[y == 1], y[y == 1]
    min_samples = min(len(X_0), len(X_1))
    # 直接截断多数类样本到少数类数量
    X_0_balanced = X_0[:min_samples] if len(X_0) > min_samples else X_0
    y_0_balanced = y_0[:min_samples] if len(y_0) > min_samples else y_0
    X_1_balanced = X_1[:min_samples] if len(X_1) > min_samples else X_1
    y_1_balanced = y_1[:min_samples] if len(y_1) > min_samples else y_1
    # 合并并打乱数据
    X_balanced = np.vstack([X_0_balanced, X_1_balanced])
    y_balanced = np.hstack([y_0_balanced, y_1_balanced])
    shuffle_idx = np.random.RandomState(random_state).permutation(len(X_balanced))
    return X_balanced[shuffle_idx], y_balanced[shuffle_idx]


# ===================== 7. 风险量化函数 =====================
def cosine_distance(v1, v2):
    """计算余弦距离：1 - 余弦相似度"""
    v1_norm = v1 / np.linalg.norm(v1, axis=1, keepdims=True)
    v2_norm = v2 / np.linalg.norm(v2, axis=1, keepdims=True)
    cos_sim = np.dot(v1_norm, v2_norm.T)
    return 1 - cos_sim


def calculate_risk_score(model, X_train_balanced, y_train_balanced, X_test, device="cpu"):
    """计算测试集样本的风险得分"""
    model.eval()
    with torch.no_grad():
        # 训练集全连接层向量
        train_embeds = torch.tensor(X_train_balanced, dtype=torch.float32).to(device)
        model(train_embeds)
        train_fc_vecs = model.fc_features.cpu().numpy()
        # 测试集全连接层向量
        test_embeds = torch.tensor(X_test, dtype=torch.float32).to(device)
        model(test_embeds)
        test_fc_vecs = model.fc_features.cpu().numpy()
    # 分离训练集人/禽向量
    train_fc_0 = train_fc_vecs[y_train_balanced == 0]
    train_fc_1 = train_fc_vecs[y_train_balanced == 1]
    # 计算每个测试样本的风险得分
    risk_scores = []
    for test_vec in test_fc_vecs:
        dist_0 = cosine_distance(test_vec.reshape(1, -1), train_fc_0).min()
        dist_1 = cosine_distance(test_vec.reshape(1, -1), train_fc_1).min()
        score_0 = 1 - (dist_0 / (dist_0 + dist_1))
        score_1 = 1 - (dist_1 / (dist_0 + dist_1))
        risk_scores.append(softmax(torch.tensor([score_0, score_1]), dim=0).numpy())
    return np.array(risk_scores)


# ===================== 8. 绘制混淆矩阵和ROC曲线 =====================
#def plot_confusion_matrix_and_roc(y_true, y_pred, y_score, save_path="vBERT_based_human_avian_classifier_confusion_matrix_roc_curve.png"):
def plot_confusion_matrix_and_roc(y_true, y_pred, y_score, save_path="vBERT_based_human_avian_classifier_confusion_matrix_roc_curve.svg"):
    """
    绘制混淆矩阵和ROC曲线的组合图
    y_true: 真实标签
    y_pred: 预测标签
    y_score: 预测概率（正类）
    """
    # 设置绘图风格
    plt.style.use('default')

    # 创建2x1的子图布局
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

    # 1. 绘制混淆矩阵
    cm = confusion_matrix(y_true, y_pred)
    cm_normalized = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]

    # 蓝白色调色板
    cmap = plt.cm.Blues

    # 绘制热图
    im = ax1.imshow(cm_normalized, interpolation='nearest', cmap=cmap)
    ax1.set_title('Confusion Matrix', fontsize=20, fontweight='bold', pad=20)

    # 设置类别标签
    classes = ['Avian (0)', 'Human (1)']
    tick_marks = np.arange(len(classes))
    ax1.set_xticks(tick_marks)
    ax1.set_yticks(tick_marks)
    ax1.set_xticklabels(classes, fontsize=16)
    ax1.set_yticklabels(classes, fontsize=16,rotation=90)

    # 设置轴标签
    ax1.set_xlabel('Predicted Label', fontsize=18, labelpad=15)
    ax1.set_ylabel('True Label', fontsize=18, labelpad=15)

    # 在每个格子中添加数值
    thresh = cm_normalized.max() / 2.
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            # 显示样本数量和归一化比例
            text = ax1.text(j, i, f'{cm[i, j]}\n({cm_normalized[i, j]:.2f})',
                            ha="center", va="center", color="white" if cm_normalized[i, j] > thresh else "black",
                            fontsize=16, fontweight='bold')

    # 2. 绘制ROC曲线
    fpr, tpr, _ = roc_curve(y_true, y_score[:, 1])
    roc_auc = auc(fpr, tpr)

    # 绘制ROC曲线
    ax2.plot(fpr, tpr, color='darkorange', lw=3, label=f'ROC curve (AUC = {roc_auc:.3f})')
    # 绘制随机猜测基准线
    ax2.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--', label='Random Guess')

    # 设置标题和标签
    ax2.set_title('Receiver Operating Characteristic (ROC) Curve', fontsize=20, fontweight='bold', pad=20)
    ax2.set_xlabel('False Positive Rate', fontsize=18, labelpad=15)
    ax2.set_ylabel('True Positive Rate', fontsize=18, labelpad=15)

    # 设置图例
    ax2.legend(loc="lower right", fontsize=14)

    # 设置轴刻度
    ax2.tick_params(axis='both', which='major', labelsize=14)

    # 调整布局
    plt.tight_layout()

    # 保存图片
    plt.savefig(save_path, dpi=600, bbox_inches='tight')
    plt.close()

    print(f"混淆矩阵和ROC曲线已保存至 {save_path}")

    # 打印分类报告
    print("\n分类报告:")
    print(classification_report(y_true, y_pred, target_names=['Avian', 'Human']))
    print(f"AUC值: {roc_auc:.3f}")

    return roc_auc


# ===================== 9. 主流程 =====================
if __name__ == "__main__":

    SEED = 42
    np.random.seed(SEED)
    torch.manual_seed(SEED)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(SEED)
        torch.cuda.manual_seed_all(SEED)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

    # 配置参数（按需修改）
    CSV_PATH = "H13579_PlosOne_aligned_and_sampled_sequences_without_empty_98_with_origin_seqid_and_seq.csv"  # 仅含seq_id的CSV路径
    NPY_PATH = "H13579_PlosOne_aligned_and_sampled_sequences_without_empty_98_with_origin_seqid_and_seq.npy"  # 单独存储的embedding npy文件
    SEQ_ID_COL = "seq_id"  # seq_id列名
    BATCH_SIZE = 64  # 增大batch提升训练速度
    EPOCHS = 10  # 训练轮数
    LR = 0.01  # 学习率
    DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

    # 数据预处理
    X_train, y_train, X_test, y_test, test_seq_ids, test_seqs = load_and_preprocess_data(
        CSV_PATH, NPY_PATH, SEQ_ID_COL, random_state=SEED
    )
    # 下采样平衡训练集
    X_train_balanced, y_train_balanced = downsample_balance(X_train, y_train, random_state=SEED)
    # 特征标准化（线性模型必需）
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train_balanced)
    X_test_scaled = scaler.transform(X_test)

    # 构建数据加载器
    train_dataset = SeqDataset(X_train_scaled, y_train_balanced)
    train_loader = DataLoader(train_dataset, batch_size=BATCH_SIZE, shuffle=True)

    # 初始化模型
    embedding_dim = X_train_scaled.shape[1]
    model = UltraSimpleLinearClassifier(embedding_dim, dropout_rate=0.1).to(DEVICE)
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.SGD(model.parameters(), lr=LR)  # SGD提速

    # 训练模型并记录loss
    trained_model, _ = fast_train_with_loss_record(
        model, train_loader, criterion, optimizer, EPOCHS, DEVICE
    )

    # 保存模型
    torch.save({
        "model_state_dict": trained_model.state_dict(),
        "scaler": scaler,
        "embedding_dim": embedding_dim
    }, "fast_linear_classifier.pth")
    print("模型已保存至 fast_linear_classifier.pth")

    # 预测与风险计算
    risk_scores = calculate_risk_score(
        trained_model, X_train_scaled, y_train_balanced, X_test_scaled, DEVICE
    )
    y_pred = np.argmax(risk_scores, axis=1)

    # 保存预测结果
    result_df = pd.DataFrame({
        "seq_id": test_seq_ids,
        "seq": test_seqs,
        "pred_label": y_pred,
        "true_label": y_test,
        "risk_score_avian": risk_scores[:, 0],
        "risk_score_human": risk_scores[:, 1]
    })
    result_df.to_csv("Predicted_prob_H13579_PlosOne_aligned_cdhit98_test_seq_split_random30.csv", index=False)
    print("预测结果已保存至 Predicted_prob_H13579_PlosOne_aligned_cdhit98_test_seq_split_random30.csv")

    # 绘制混淆矩阵和ROC曲线
    plot_confusion_matrix_and_roc(y_test, y_pred, risk_scores)
