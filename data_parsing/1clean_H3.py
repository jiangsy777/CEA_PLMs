from Bio import SeqIO

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


def process_fasta(input_file, output_file, ref_name):
    fasta_file = parse_fasta(input_file)
    seq_id_lst = fasta_file[0]
    seq_seq_lst = fasta_file[1]
    # 读取FASTA文件，存储序列名和序列
    seq_dict = dict(zip(seq_id_lst,seq_seq_lst))
    print(list(seq_dict.keys())[1])
    ref_seq = seq_dict[ref_name]
    seq_len = len(ref_seq)
    all_seqs = list(seq_dict.values())
    seq_names = list(seq_dict.keys())
    print(len(seq_names))

    # 标记参考序列中的gap列（空格位置）
    gap_cols = [i for i, char in enumerate(ref_seq) if char == "-"]
    keep_cols = []  # 最终保留的列索引
    delete_cols = []

    # 计算每个gap列的空格比例，筛选出需保留的列
    for col in gap_cols:
        gap_count = sum(1 for seq in all_seqs if seq[col] == "-")
        gap_ratio = gap_count / len(all_seqs)
        if gap_ratio < 0.95:  # 仅保留空格比例<95%的列
            keep_cols.append(col)
        else:
            delete_cols.append(col)
    print(keep_cols)
    print(delete_cols)

    # 保留所有非gap列 + 筛选后的gap列
    all_cols = set(range(seq_len))
    final_keep_cols = sorted((all_cols - set(gap_cols)) | set(keep_cols))

    # 过滤序列
    filtered = []
    for name, seq in zip(seq_names, all_seqs):
        # 检查该序列在“待删列”中是否有非空格（有则删除）
        has_non_gap_in_delete_cols = any(seq[col] != "-" for col in delete_cols)
        if not has_non_gap_in_delete_cols:
            # 仅对符合条件的序列，提取“保留列”的内容
            trimmed = "".join([seq[col] for col in final_keep_cols])
            filtered.append((name, trimmed))
    n = 0
    # 写入输出文件
    with open(output_file, "w") as f:
        ref_name1 = 'A_Aichi_2_1968_H3N2_align'
        ref_seq1 = 'MKTIIALSYIFCLALGQDLPGNDNSTATLCLGHHAVPNGTLVKTITDDQIEVTNATELVQSSSTGKICNN-PHRILDGIDCTLIDALLGDPHCDVFQNE-TWDLFVERSK-AFSNCYPYDVPDYASLRSLVASSGTLEFITEGF---TWTGVTQN-GGSNACKR-GPGSGFFSRLNWLTKSG--STYPVLNVTMPNNDNFDKLYIWGIHHPSTNQEQTSLYVQASGRVTVSTRRSQQTIIPNIGSRPWVRGLSSRISIYWTIVKPGDVLVINSNGNLIAPRGYFKMRT----------GKSSIMRSDAPIDTCISECITPNGSIPNDKPFQNVNKITYGACPKYVKQNTLKLATGMRNVPEKQT----RGLFGAIAGFIENGWEGMIDGWYGFRHQNSEGTGQAADLKSTQAAIDQINGKLNRVIEKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTRRQLRENAEDMGNGCFKIYHKCDNACIESIRNGTYDHDVYRDEALNNRFQIKGVELKSGYKDWILWISFAISCFLLCVVLLGFIMWACQRGNIRCNICI'
        ref_seq1 = ref_seq1.replace(' ','')
        f.write(f">{ref_name1}\n{ref_seq1}\n")

        ref_name2 = 'A_Aichi_2_1968_H3N2'
        ref_seq2 = ref_seq1.replace('-','')
        f.write(f">{ref_name2}\n{ref_seq2}\n")

        for name, seq in filtered:
            n += 1
            f.write(f">{name}\n{seq}\n")

    print(n)
if __name__ == "__main__":
    # 替换为你的输入输出文件路径



    process_fasta("../data/mafft_aligned/HA_with_aa_seq_without_wrong_seq_H3_mafft.fasta", "../data/generate/HA_with_aa_seq_without_wrong_seq_H3_mafft_del_empty_with_plosone.fasta",ref_name='A/England/220480330/2022_Human_H3')





