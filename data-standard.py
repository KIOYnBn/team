# 本文件用于将原始数据整合成本模型使用的标准数据格式

import os
from Bio import SeqIO

# 定义输入输出文件路径
input_file = "input.fasta"
output_file = "output.fasta"
cd_hit_path = "/path/to/cd-hit"

# 运行CD-HIT命令
os.system(f"{cd_hit_path} -i {input_file} -o {output_file} -c 0.9 -n 5")

# 读取去冗余后的序列
sequences = SeqIO.parse(output_file, "fasta")
for seq_record in sequences:
    print(f">{seq_record.id}\n{seq_record.seq}")


