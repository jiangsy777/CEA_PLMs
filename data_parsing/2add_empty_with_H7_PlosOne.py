def parse_fasta(input_file):
    """自定义FASTA解析，返回 [(序列名, 完整序列), ...]"""
    fasta_list = []
    current_name, current_seq = "", []
    with open(input_file, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name:
                    fasta_list.append((current_name, "".join(current_seq)))
                current_name = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_name:
            fasta_list.append((current_name, "".join(current_seq)))
    return fasta_list


def write_fasta(fasta_list, output_file):
    """自定义FASTA写入，适配常规格式"""
    with open(output_file, "w", encoding="utf-8") as f:
        for name, seq in fasta_list:
            f.write(f">{name}\n{seq}\n")


def align_to_first_seq(input_fasta, output_fasta, start_loc, end_loc, gap_char="-"):

    # 1. 解析FASTA，分离模板、基准及目标序列
    fasta_list = parse_fasta(input_fasta)
    if len(fasta_list) < 2:
        raise ValueError("FASTA文件序列不足2条，无法执行对齐")

    template_name, template_seq = fasta_list[0]
    base_name, base_seq = fasta_list[1]
    target_seqs = fasta_list[1:]  # 第二条及以后需插空格

    # 2. 定位基准序列原有空格，筛选模板需插入的空格索引
    template_seq_list = list(template_seq)
    base_seq_list = list(base_seq)
    # base_gap_set = set(idx for idx, c in enumerate(base_seq_list) if c == gap_char)
    template_insert_indices0 = [
        idx for idx, c in enumerate(template_seq_list)
        if c == gap_char]
    template_insert_indices = []
    processed_list = [(template_name, template_seq[start_loc:end_loc])]
    for insert_idx in template_insert_indices0:
        if base_seq_list[insert_idx] != gap_char:
            base_seq_list.insert(insert_idx, gap_char)
            template_insert_indices.append(insert_idx)
    base_seq_list1 = base_seq_list[start_loc:end_loc]
    processed_list.append((base_name, "".join(base_seq_list1)))


    # 3. 从前往后给目标序列插入空格

    target_seqs1 = target_seqs[1:]
    for tar_name, tar_seq in target_seqs1:
        tar_seq_list = list(tar_seq)
        for insert_idx in template_insert_indices:
            tar_seq_list.insert(insert_idx, gap_char)
        tar_seq_list1 = tar_seq_list[start_loc:end_loc]
        processed_list.append((tar_name, "".join(tar_seq_list1)))

    # 4. 写入结果
    write_fasta(processed_list, output_fasta)
    print(f"对齐完成！共处理{len(target_seqs)}条序列，结果保存至{output_fasta}")


if __name__ == "__main__":
    # 替换为你的文件路径，无需安装任何依赖，直接运行
    INPUT = "../data/generate/HA_with_aa_seq_without_wrong_seq_H7_mafft_del_empty_with_plosone.fasta"
    OUTPUT = "../data/generate/HA_with_aa_seq_without_wrong_seq_H7_mafft_del_empty_with_plosone_final.fasta"
    start_loc1 = 19
    end_loc1 = 361
    start_loc = start_loc1-1
    end_loc = end_loc1
    align_to_first_seq(input_fasta=INPUT, output_fasta=OUTPUT,start_loc=start_loc,end_loc=end_loc)


