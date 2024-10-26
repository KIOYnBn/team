from Bio.PDB import *
import os
from icecream import ic

# Bio初始化
pdbl = PDBList()

# 文件夹创建
os.makedirs('./pdb_sequence', exist_ok=True)

## 数据下载
pdb_lst = ["3JCL", "2FAT"] # 目标列表

# 开始获取
for count in range(len(pdb_lst)):
    pdbl.retrieve_pdb_file(pdb_lst[count],pdir='.', file_format='pdb')
    file_path = f'pdb{pdb_lst[count].lower()}.ent'
    file_path_mid = f'./pdb_sequence/{pdb_lst[count]}'
    with open (file_path, "r") as f:
        with open (file_path_mid, "w+") as f_mid:
            context = f.readlines()
            for index in range(len(context)):
                judge = 'SEQRES' in context[index]
                if judge:
                    f_mid.write(context[index]+'\n')
    os.remove(file_path)
    with open(file_path_mid, "r") as f:
        context_final = f.readlines()
        protein_sequence = ""
        for count_two in range(0,len(context_final), 2):
            seq = context_final[count_two].split(" ")[7:-1]
            while "" in seq:
                seq.remove('')
            for count_three in range(len(seq)):
                protein_sequence += seq[count_three]
    with open(f"{file_path_mid}_final", "w") as f:
        f.write(protein_sequence)
    os.remove(file_path_mid)


'''
def get_sequence(structure):
    model = structure[0]  # 通常只有一个模型
    chain = model['A']  # 假设只有一个链
    sequence = ''
    for residue in chain:
        if residue.get_resname() in three_to_one:
            sequence += three_to_one[residue.get_resname()]
    return sequence


# 定义一个函数来生成接触图
def generate_contact_map(structure, threshold=5.0):
    model = structure[0]
    chain = model['A']
    contact_map = {}
    residues = list(chain.get_residues())
    for i in range(len(residues)):
        for j in range(i + 1, len(residues)):
            # 计算两个残基之间的最小原子距离
            min_dist = float('inf')
            for atom1 in residues[i]:
                for atom2 in residues[j]:
                    dist = atom1 - atom2
                    min_dist = min(min_dist, dist)
            if min_dist < threshold:
                contact_map[(residues[i].get_id()[1], residues[j].get_id()[1])] = min_dist
    return contact_map



    # 解析结构
    structure = parser.get_structure(pdb_id, file_path)

    

    # 生成接触图
    contact_map = generate_contact_map(structure)
    print(f"Contact map of {pdb_id}: {contact_map}")

'''
