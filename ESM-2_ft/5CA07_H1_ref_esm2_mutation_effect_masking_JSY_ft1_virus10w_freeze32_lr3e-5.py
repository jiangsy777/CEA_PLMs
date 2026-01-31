from transformers import AutoModelForMaskedLM, AutoTokenizer
import torch
import torch.nn.functional as F
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os


def esm2_direct_logit_comparison(sequence, position, mutant_aa, model, tokenizer, device):
    """
    直接比较野生型和突变型的logits来计算突变效应
    这种方法更精确，避免掩码预测的偏差
    """
    # 准备野生型序列
    inputs = tokenizer(sequence, return_tensors="pt", padding=True, truncation=True)
    inputs = {key: value.to(device) for key, value in inputs.items()}

    # 获取所有位置的logits
    with torch.no_grad():
        outputs = model(**inputs)

    logits0 = outputs.logits[0]  # [序列长度, 词汇表大小33]
    logits = F.log_softmax(logits0, dim=-1)  # [序列长度, 词汇表大小33]

    #print(logits.shape)

    # 调整位置索引 (考虑特殊的起始token)
    adjusted_position = position + 1  # +1 因为第一个token是<cls>

    wt_aa = sequence[position]
    wt_token_id = tokenizer.convert_tokens_to_ids(wt_aa)  # 野生的token ID
    mutant_token_id = tokenizer.convert_tokens_to_ids(mutant_aa)  # 突变的token ID

    # 验证token ID是否有效
    if wt_token_id == tokenizer.unk_token_id:
        raise ValueError(f"无效的野生型氨基酸: {wt_aa}")
    if mutant_token_id == tokenizer.unk_token_id:
        raise ValueError(f"无效的突变型氨基酸: {mutant_aa}")

    # 获取野生型和突变型的logits
    wt_logit = logits[adjusted_position, wt_token_id].item()
    mutant_logit = logits[adjusted_position, mutant_token_id].item()

    # 计算对数概率比
    log_prob_ratio = mutant_logit - wt_logit

    return {
        "wildtype_aa": wt_aa,
        "mutant_aa": mutant_aa,
        "position": position,
        "log_prob_ratio": log_prob_ratio,
        "wt_logit": wt_logit,
        "mutant_logit": mutant_logit
    }


def analyze_all_mutations_for_sequence(sequence, model, tokenizer, device, output_csv="mutation_effects.csv"):
    """
    分析参考序列每个位点的氨基酸突变为20种氨基酸的效应
    将结果保存到CSV文件
    """
    # 20种标准氨基酸
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"

    results = []

    print(f"开始分析序列，长度: {len(sequence)}")
    print(f"设备: {device}")

    for position in range(len(sequence)):
        wt_aa = sequence[position]
        print(f"处理位置 {position + 1}/{len(sequence)}: {wt_aa}")

        # 对每个位置，突变为20种氨基酸
        for mutant_aa in amino_acids:
            # if mutant_aa == wt_aa:
            #     continue  # 跳过相同的氨基酸

            try:
                result = esm2_direct_logit_comparison(sequence, position, mutant_aa, model, tokenizer, device)

                # 记录结果
                mutation_id = f"{wt_aa}{position + 1}{mutant_aa}"  # 1-based位置用于显示
                results.append({
                    "position": position,  # 0-based位置
                    "mutation": mutation_id,
                    "wt_aa": wt_aa,
                    "mutant_aa": mutant_aa,
                    "delta_logit": result["log_prob_ratio"],
                    "wt_logit": result["wt_logit"],
                    "mutant_logit": result["mutant_logit"]
                })

            except Exception as e:
                print(f"处理突变 {position}{wt_aa}→{mutant_aa} 时出错: {e}")
                continue

    # 创建DataFrame并保存为CSV
    df = pd.DataFrame(results)

    # 按位置排序
    #df = df.sort_values("position")

    # 保存到CSV
    df.to_csv(output_csv, index=False)
    print(f"\n结果已保存到: {output_csv}")
    print(f"总共分析了 {len(df)} 个突变")

    return df


def load_model_and_tokenizer_local(model_path):
    """
    从本地路径加载模型和tokenizer
    """
    print(f"从本地路径加载模型: {model_path}")

    # 检查模型路径是否存在
    if not os.path.exists(model_path):
        raise FileNotFoundError(f"模型路径不存在: {model_path}")

    # 加载tokenizer
    tokenizer = AutoTokenizer.from_pretrained(model_path, local_files_only=True)
    vocab = tokenizer.get_vocab()
    print(vocab)

    # 加载模型
    model = AutoModelForMaskedLM.from_pretrained(model_path, local_files_only=True)

    return model, tokenizer


# 使用示例
if __name__ == "__main__":
    # 设置设备
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"使用设备: {device}")
    #checkpoints = [i for i in [1000,1180]]
    checkpoints = [1800,2000,2200]

    for cp in checkpoints:
        print('cp=',cp)
        cp1 = str(cp)
    
        # 本地模型路径 - 请根据实际情况修改
        local_model_path = "./esm2_finetuned_virus10w_freeze32layers_lr3e-5/checkpoint-"+cp1   # 替换为你的实际路径

        # 加载模型和tokenizer
        try:
            model, tokenizer = load_model_and_tokenizer_local(local_model_path)
            model.to(device)
            #tokenizer.to(device)
            model.eval()
            print("模型加载成功!")
        except Exception as e:
            print(f"模型加载失败: {e}")
            exit(1)

        # 示例野生蛋白质序列 (TP53片段)
        #tp53_sequence = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD"
        #tp53_sequence = "MEEPQSD"
        #tp53_sequence = 'KAKLLVLLYAFVATDADTICIGYHANNSTDTVDTILEKNVAVTHSVNLLEDSHNGKLCKLKGIAPLQLGKCNITGWLLGNPECDSLLPARSWSYIVETPNSENGACYPGDLIDYEELREQLSSVSSLERFEIFPKESSWPNHTFNGVTVSCSHRGKSSFYRNLLWLTKKGDSYPKLTNSYVNNKGKEVLVLWGVHHPSSSDEQQSLYSNGNAYVSVASSNYNRRFTPEIAARPKVRDQHGRMNYYWTLLEPGDTIIFEATGNLIAPWYAFALSRGFESGIITSNASMHECNTKCQTPQGAINSNLPFQNIHPVTIGECPKYVRSTKLRMVTGLRNIPSIQYRGLFGAIAGFIEGGWTGMIDGWYGYHHQNEQGSGYAADQKSTQNAINGITNKVNSVIEKMNTQFTAVGKEFNNLEKRMENLNKKVDDGFLDIWTYNAELLVLLENERTLDFHDLNVKNLYEKVKSQLKNNAKEIGNGCFEFYHKCDNECMESVRNGTYDYPKYSEESKLNREKIDGVKLESMGVYQILAIYSTVASSLVLLVSLGAISFWMCSNGSLQCRICI'
        #tp53_sequence = 'NITNLCPFDEVFNATRFASVYAWNRKRISNCVADYSVLYNLAPFFTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGNIADYNYKLPDDFTGCVIAWNSNKLDSKVSGNYNYLYRLFRKSNLKPFERDISTEIYQAGNKPCNGVAGFNCYFPLRSYSFRPTYGVGHQPYRVVVLSFELLHAPATVCGPKKST'
        tp53_sequence = 'MKAILVVLLYTFATANADTLCIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDKHNGKLCKLRGVAPLHLGKCNIAGWILGNPECESLSTASSWSYIVETPSSDNGTCYPGDFIDYEELREQLSSVSSFERFEIFPKTSSWPNHDSNKGVTAACPHAGAKSFYKNLIWLVKKGNSYPKLSKSYINDKGKEVLVLWGIHHPPTSADQQSLYQNADAYVFVGSSRYSKKFKPEIAIRPKVRDQEGRMNYYWTLVEPGDKITFEATGNLVVPRYAFAMERNAGSGIIISDTPVHDCNTTCQTPKGAINTSLPFQNIHPITIGKCPKYVKSTKLRLATGLRNIPSIQSRGLFGAIAGFIEGGWTGMVDGWYGYHHQNEQGSGYAADLKSTQNAIDEITNKVNSVIEKMNTQFTAVGKEFNHLEKRIENLNKKVDDGFLDIWTYNAELLVLLENERTLDYHDSNVKNLYEKVRSQLKNNAKEIGNGCFEFYHKCDNTCMESVKNGTYDYPKYSEEAKLNREEIDGVKLESTRIYQILAIYSTVASSLVLVVSLGAISFWMCSNGSLQCRICI'

        print("=== ESM2突变效应分析 (本地模型版本) ===")
        print(f"序列长度: {len(tp53_sequence)}")

        # 分析所有可能的突变
        results_df = analyze_all_mutations_for_sequence(
            sequence=tp53_sequence,
            model=model,
            tokenizer=tokenizer,
            device=device,
            output_csv="./CA07_ref_HA_logist_new_virus10w_freeze32_lr3e-5/HA_mutation_effects_esm2_ft"+cp1+".csv"
        )

        # 显示前几行结果
        print("\n前10个突变结果:")
        print(results_df.head(10))

        # 统计信息
        print(f"\n统计信息:")
        print(f"总突变数: {len(results_df)}")
        print(f"Delta Logit 范围: [{results_df['delta_logit'].min():.3f}, {results_df['delta_logit'].max():.3f}]")
        print(f"平均 Delta Logit: {results_df['delta_logit'].mean():.3f}")

        # 有害突变统计 (假设delta_logit < -2为有害)
        deleterious_mutations = results_df[results_df['delta_logit'] < -2]
        print(f"有害突变数 (delta_logit < -2): {len(deleterious_mutations)}")
