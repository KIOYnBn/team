from icecream import ic
import os

def seq_integrate(uniprot_id):
    os.makedirs("./sequence_uni_all", exist_ok=True)
    with open("./sequence_uni_all/all_sequence.fasta", "a+") as all_file:
            with open('./sequence_uni' + "/" + f'{uniprot_id}.fasta', 'r') as pre_file:
                all_file.write(pre_file.read())


with open("uniprot_pdb", "r") as f:
    con = f.readlines()
ic(con)

for count in range(1, 10):
    a = con[count].split(",")[0]
    ic(a)
    seq_integrate(a)
