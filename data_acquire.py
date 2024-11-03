from Bio.PDB import *
import os
from pyuniprot import Uniprot
from black.comments import children_contains_fmt_on
from icecream import ic
from rdkit import Chem
from collections import defaultdict
import numpy as np
from time import sleep, strftime, localtime
import requests
from Bio import SeqIO
from io import StringIO
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import pairwise_distances


def url_pdb_fasta(pdbid):
    url = f"https://www.rcsb.org/fasta/entry/{pdbid}/display"
    # 发送HTTP获取FASTA数据
    response = requests.get(url)
    try:
        with open(f"./{pdbid}", "w") as f:
            fasta_records = list(SeqIO.parse(StringIO(response.text), "fasta"))
            for record in fasta_records:
                print(f">{record.id} {record.description}")
                print(record.seq)
                f.write(f">{record.id} {record.description}\n{str(record.seq)}")
    except Exception as e:
        with open(f"./error_log/url_pdb_fasta{pdbid}_error", "a") as f:
            f.write(f"error: {e} to  {pdbid}")


def seq_pdb(pdb_lst, target_name, sequence_judge = True):
    # Bio初始化
    pdbl = PDBList()
    parser = PDBParser(QUIET=True)
    ppb = PPBuilder()
    # 文件夹创建
    #os.makedirs('./pdb_sequence', exist_ok=True)
    for count in range(len(pdb_lst)):
        pdbl.retrieve_pdb_file(pdb_lst[count],pdir='.', file_format='pdb', overwrite=True)
        file_path = f'pdb{pdb_lst[count].lower()}.ent'
        try:
            if sequence_judge:# 隐性序列获取（fasta）
                structure = parser.get_structure(file=file_path, id=None)
                chains = structure.get_chains()
                chain_ids = [chain.get_id() for chain in chains]  # 链的名称
                chain_id = chain_ids[0] 
                contact_map = protein_contact_map(pdb=pdb_lst[count], chain=chain_id) # 
                
                sequence = str(ppb.build_peptides(structure)[0].get_sequence())
                ic(sequence)
                '''# 全部链的信息
                for count_pp in range(len(pp_all)):
                    sequence = pp_all[count_pp].get_sequence()
                    ic(sequence)
                    all_sequence += ' '+str(sequence) + "\n"
                '''
                # with open(f"./DUD-E/all/{target_name[count]}/{pdb_lst[count]}", "w") as f:
                  #  f.write(">"+pdb_lst[count]+"\n"+sequence+"\n")
                np.save(f"./DUD-E/all/{target_name[count]}/{pdb_lst[count]}_contact_map", contact_map)
                os.remove(file_path)
                sleep(0.1)

            else:# 显性序列获取（seqres）
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
        except Exception as e:
            with open("error_log/acquire_seq_pdb_err.log.txt", "a") as f:
                f.write(f"{time_now}: an error occured for {pdb_lst[count]}: {e}\n")
            continue


def protein_contact_map(pdb, chain, seq_gap=4, contact_cutoff=6, map=False): # pdb和chain输入均为字符串
    # map是指是否画图
    # 若没有pdb本地文件, 使用这个
    '''
    url = f'https://files.rcsb.org/download/{pdb}.pdb'
    r = requests.get(url)
    r.raise_for_status()
    lines = r.text.splitlines()
    '''
    with open(f"./pdb{pdb}.ent", "r") as f: # 若有本地pdb文件, 使用这个
        lines = f.readlines()
    out = []
    for line in lines:
        if line.startswith('ATOM ') and line.split()[4] == chain and line.split()[2] == 'CA':
            resi = line.split()[5]
            resn = line.split()[3]
            x = line.split()[6]
            y = line.split()[7]
            z = line.split()[8]
            out.append([resi, resn, x, y, z])
    df = pd.DataFrame(out, columns=['res_num', 'res_name', 'x', 'y', 'z'])
    dist_arr = pairwise_distances(df[['x', 'y', 'z']].values)
    if map:
        fig = plt.figure(figsize=(6, 5))
        p = plt.imshow(dist_arr, cmap='viridis_r', vmax=24, vmin=0)
        plt.colorbar(p, label='Distance ($\mathrm{\AA}$)')
        plt.xlabel('Residue i')
        plt.ylabel('Residue j')
        plt.title('Luciferase (1LUC) Distrogram')
        plt.savefig(f'{pdb}_distro.png', bbox_inches='tight', dpi=300)
    return dist_arr


def seq_uniprot(uniprot_id, protein_name="example"):
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
            f.write(f"{time_now}: an error occured for {uniprot_id}: {e}\n")


def seq_integrate(uniprot_id):
    os.makedirs("./sequence_uni_all", exist_ok=True)
    try:
        with open("./sequence_uni_all/all_sequence.fasta", "a+") as all_file:
                with open('./sequence_uni' + "/" + f'{uniprot_id}.fasta', 'r') as pre_file:
                    all_file.write(pre_file.read())
    except Exception as e:
        with open("error_log/acquire_seq_integrate_err.log.txt", "a") as f:
            f.write(f"an error occured all_seq for {uniprot_id}: {e}\n")


def create_atoms(mol):
    atom_dict = defaultdict(lambda: len(atom_dict))
    """Create a list of atom (e.g., hydrogen and oxygen) IDs
    considering the aromaticity."""
    atoms = [a.GetSymbol() for a in mol.GetAtoms()]
    for a in mol.GetAromaticAtoms():
        i = a.GetIdx()
        atoms[i] = (atoms[i], 'aromatic')
    atoms = [atom_dict[a] for a in atoms]
    return np.array(atoms)


def create_ijbonddict(mol):
    """Create a dictionary, which each key is a node ID
    and each value is the tuples of its neighboring node
    and bond (e.g., single and double) IDs."""
    bond_dict = defaultdict(lambda: len(bond_dict))
    i_jbond_dict = defaultdict(lambda: [])
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        bond = bond_dict[str(b.GetBondType())]
        i_jbond_dict[i].append((j, bond))
        i_jbond_dict[j].append((i, bond))
    return i_jbond_dict


def extract_fingerprints(atoms, i_jbond_dict, radius):
    """Extract the r-radius subgraphs (i.e., fingerprints)
    from a molecular graph using Weisfeiler-Lehman algorithm."""
    fingerprint_dict = defaultdict(lambda: len(fingerprint_dict))
    edge_dict = defaultdict(lambda: len(edge_dict))
    if (len(atoms) == 1) or (radius == 0):
        fingerprints = [fingerprint_dict[a] for a in atoms]
    else:
        nodes = atoms
        i_jedge_dict = i_jbond_dict
        for _ in range(radius):
            """Update each node ID considering its neighboring nodes and edges
            (i.e., r-radius subgraphs or fingerprints)."""
            fingerprints = []
            for i, j_edge in i_jedge_dict.items():
                neighbors = [(nodes[j], edge) for j, edge in j_edge]
                fingerprint = (nodes[i], tuple(sorted(neighbors)))
                fingerprints.append(fingerprint_dict[fingerprint])
            nodes = fingerprints
            """Also update each edge ID considering two nodes
            on its both sides."""
            _i_jedge_dict = defaultdict(lambda: [])
            for i, j_edge in i_jedge_dict.items():
                for j, edge in j_edge:
                    both_side = tuple(sorted((nodes[i], nodes[j])))
                    edge = edge_dict[(both_side, edge)]
                    _i_jedge_dict[i].append((j, edge))
            i_jedge_dict = _i_jedge_dict

    return np.array(fingerprints)


def create_adjacency(mol):
    adjacency = Chem.GetAdjacencyMatrix(mol)
    return np.array(adjacency)


time_now = strftime("%Y-%m-%d %H:%M:%S", localtime())
