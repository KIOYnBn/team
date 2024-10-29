import subprocess

# 定义cd-hit的输入文件和输出文件
input_file = "sequences.fasta"
output_file = "sequences_cdhit.fasta"
cdhit_db = "cdhit_database"

# 定义cd-hit的参数
cdhit_options = {
    '-i': input_file,         # 输入文件
    '-o': output_file,        # 输出文件
    '-c': '0.9',              # 序列相似性阈值
    '-M': '16000',            # 内存限制
    '-T': '4',                # 线程数
    '-d': '0'                 # 0表示核酸，1表示蛋白质
}

# 将参数转换为命令行参数列表
cdhit_command = ["cd-hit"]
for option, value in cdhit_options.items():
    cdhit_command.extend([option, value])

# 调用cd-hit
try:
    subprocess.run(cdhit_command, check=True)
    print(f"CD-HIT去冗余完成，结果保存在 {output_file}")
except subprocess.CalledProcessError as e:
    print(f"CD-HIT运行出错: {e}")
