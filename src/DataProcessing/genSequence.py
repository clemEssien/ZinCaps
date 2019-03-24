import os

from glob import glob


def capitalise(str):
    return str[0:6].upper()+str[6:len(str)]

DATASET_DIR = "../../data/dataset/"
RES_DIR = "../../data/res/"
positive_pids = []
negative_pids = []
for files in glob(DATASET_DIR+"*"):
    filename = os.path.split(files)[1]
    file = open(DATASET_DIR+filename,"r")
    for f in file:
        if "positive" in filename :
            positive_pids.append(f[0:6]);
        else:
            negative_pids.append(f[0:6])

positive_pids = list(set(positive_pids))
negative_pids = list(set(negative_pids))
combined_pids = list(set(positive_pids +negative_pids))

with open(RES_DIR+"pdb_seqres.txt", "r") as file:
   fasta_list = []
   data = file.read()
   seq_list = data.split(">")

   for pdb_id in combined_pids:
       for seq in seq_list:
           if(pdb_id.strip().lower() == seq[0:6].strip().lower()):
               with open(RES_DIR+"sequence.fasta", "a") as fasta_file:
                    fasta_file.write(capitalise(seq))