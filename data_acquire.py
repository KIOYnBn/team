from Bio.PDB import *
import os
from pyuniprot import Uniprot
from black.comments import children_contains_fmt_on
from icecream import ic


# 开始获取
def seq_pdb(pdb_lst, sequence_judge = True):
    # Bio初始化
    pdbl = PDBList()
    parser = PDBParser(QUIET=True)
    ppb = PPBuilder()
    # 文件夹创建
    os.makedirs('./pdb_sequence', exist_ok=True)
    for count in range(len(pdb_lst)):
        pdbl.retrieve_pdb_file(pdb_lst[count],pdir='.', file_format='pdb')
        file_path = f'pdb{pdb_lst[count].lower()}.ent'
        file_path_mid = f'./pdb_sequence/{pdb_lst[count]}'

        if sequence_judge:# 隐性序列获取（fasta）
            # 这里对于fasta里的序列，只取第一条链
            structure = parser.get_structure(file=file_path, id=None)
            chains = structure.get_chains()
            sequence = str(ppb.build_peptides(structure)[0].get_sequence())
            '''# 全部链的信息
            for count_pp in range(len(pp_all)):
                sequence = pp_all[count_pp].get_sequence()
                ic(sequence)
                all_sequence += ' '+str(sequence) + "\n"
            '''
            with open(f"{file_path_mid}_final", "w") as f:
                f.write(">"+pdb_lst[count]+"\n"+sequence)
            os.remove(file_path)

        else:# 显性序列获取（seqres）
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
                # chain_name = context_final[0].split(" ")[3]
                # protein_sequence = chain_name+" "
                protein_sequence = ''
                chain_name = ''
                for count_two in range(0,len(context_final), 2):
                    seq_all = context_final[count_two].split(" ")[:-1]
                    while "" in seq_all:
                        seq_all.remove("")
                    ic(seq_all)
                    if chain_name == '':
                        chain_name = seq_all[2]
                        protein_sequence += chain_name+" "
                    now_chain_name = seq_all[2]
                    if chain_name == now_chain_name:
                        pass
                    else:
                        chain_name = now_chain_name
                        protein_sequence += "\n"+now_chain_name+" "
                    for fragment in seq_all[4:]:
                        protein_sequence += fragment
            with open(f"{file_path_mid}_final", "w") as f:
                f.write(protein_sequence)
            os.remove(file_path_mid)


def seq_uniprot(uniprot_id, protein_name):
    dir_seq = './sequence_uni'
    os.makedirs(dir_seq, exist_ok=True)
    try:
        uniprot = Uniprot(uniprot_id, local_download_dir=dir_seq, save_txt=True)
        sequence = str(uniprot.category_lines['SQ']).split("'")[1]
        fasta_format = f">{uniprot_id} {protein_name}\n{sequence}\n"
        with open(dir_seq + "/" + f'{uniprot_id}.fasta', 'w') as fasta_file:
            fasta_file.write(fasta_format)
        os.remove(dir_seq + "/" + f'{uniprot_id}.txt')
    except Exception as e:
        with open("error_log/acquire_seq_uni_err.log.txt", "a") as f:
            f.write(f"an error occured for {uniprot_id}: {e}\n")


def seq_integrate(uniprot_id):
    os.makedirs("./sequence_uni_all", exist_ok=True)
    try:
        with open("./sequence_uni_all/all_sequence.fasta", "a+") as all_file:
                with open('./sequence_uni' + "/" + f'{uniprot_id}.fasta', 'r') as pre_file:
                    all_file.write(pre_file.read())
    except Exception as e:
        with open("error_log/acquire_seq_integrate_err.log.txt", "a") as f:
            f.write(f"an error occured all_seq for {uniprot_id}: {e}\n")


os.makedirs("error_log", exist_ok=True)
## 数据下载
pdb_lst = ["4B3E", "2FAT"] # 目标列表
# seq_pdb(pdb_lst)

# 设置Uniprot ID和蛋白质名称
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
