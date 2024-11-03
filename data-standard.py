import os
from data_acquire import *


# 运行错误日志
os.makedirs("error_log", exist_ok=True)


# url_pdb_fasta("3eml") # 使用url获取pdbfasta序列


'''
# DUD-E数据集 使用biopython获取蛋白质信息
with open("./DUD-E/target-protein.txt", 'r') as f:
    pdb_lst = [line.strip() for line in f.readlines()] # 目标列表
    ic(pdb_lst)
target_name, protein = [], []
for count in pdb_lst:
    all = count.split(" ")
    target_name.append(all[0])
    protein.append(all[1])
ic(target_name, protein)
seq_pdb(protein[1:], target_name[1:])
'''


'''
# 获取uniprot的蛋白质
with open("uniprot_pdb", "r") as f:
    context = f.readlines()
for count in range(1, len(context)):
    uniprot_ids = context[count].split(",")[0]
    seq_uniprot(uniprot_ids, "example protein")
    seq_integrate(uniprot_ids)

uniprot_id = "P01116"
protein_name = "Example Protein"
# seq_uniprot(uniprot_ids, protein_name)
'''


# 分子部分
smiles = "C#CCOc3nc(c1ccccc1)nc4sc2CCCc2c34" # 分子的smile representation
mol = Chem.AddHs(Chem.MolFromSmiles(smiles))  # Consider hydrogens.
atoms = create_atoms(mol)
a = create_adjacency(mol)
i_jbond_dict = create_ijbonddict(mol)
fingerprints = extract_fingerprints(atoms, i_jbond_dict, radius=1) # radius即半径, 自定义
ic(atoms, i_jbond_dict, fingerprints, a)

