from Bio.PDB import *
from Bio import SeqIO

# Bio初始化
pdbl = PDBList()
parser = PDBParser(QUIET=True)
ppb = PPBuilder()
## 数据下载
pdb_lst = ["5FQD", "2FAT"] # 目标列表


for count in range(len(pdb_lst)):
    pdbl.retrieve_pdb_file(pdb_lst[count],pdir='.', file_format='pdb')
    file_path = f'pdb{pdb_lst[count].lower()}.ent'
    file_path_mid = f'./pdb_sequence/{pdb_lst[count]}'
    with open(file_path) as handle:
        sequence = next(SeqIO.parse(handle, "pdb-atom"))
    with open(file_path, "w") as output_handle:
        SeqIO.write(sequence, output_handle, "fasta")