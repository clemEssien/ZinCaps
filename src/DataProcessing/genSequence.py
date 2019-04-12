from tqdm import tqdm
import os
from glob import glob
from collections import defaultdict
import re

def capitalise(str):
    return ">" + str[0:6].upper() + str[6:len(str)]

def f_comma(my_str, group, char):
    my_str = str(my_str)
    return char.join(my_str[i:i+group] for i in range(0, len(my_str), group))

def index_of_xters_in_string(string, xter):
	return [match.start() for match in re.finditer(re.escape(xter), string)]

def insert(string, array,xter):
  sequence= [0];
  new_string = "";
  i = 1;
  for i in range(len(array)):
    sequence.append(array[i])
    new_string += string[sequence[i]:sequence[i+1]]+xter
  new_string += string[sequence[-1]:]+"\n"
  return new_string

DATASET_DIR = "../../data/dataset/"
RES_DIR = "../../data/resource/"
positive_pids = []
negative_pids = []
temp_pos_list = []
temp_neg_list = []
pos_dict = defaultdict(list)
neg_dict = defaultdict(list)

for files in glob(DATASET_DIR+"*"):
    filename = os.path.split(files)[1]
    file = open(DATASET_DIR+filename,"r")
    for f in file:
        if "positive" in filename :
            id = f[0:6]
            positive_pids.append(id);
            if positive_pids[len(positive_pids)-1] == id:
                temp_pos_list.append(f[11:14].strip())
            pos_dict[id].append(temp_pos_list)
            temp_pos_list = []
        else:
            id = f[0:6]
            negative_pids.append(id);
            if negative_pids[len(negative_pids)-1] == id:
                temp_neg_list.append(f[11:14].strip())
            neg_dict[id].append(temp_neg_list)
            temp_neg_list = []

positive_pids = list(set(positive_pids))
negative_pids = list(set(negative_pids))
combined_pids = list(set(positive_pids +negative_pids))

print(pos_dict['3EB5_A'])
print(type(pos_dict['3EB5_A']))
print(neg_dict)
with open(RES_DIR+"pdb_seqres.txt", "r") as file:
   fasta_list = []
   data = file.read()
   seq_list = data.split(">")

   for pdb_id in tqdm(positive_pids):
       for seq in tqdm(seq_list):
           if pdb_id.strip().lower() == seq[0:6].strip().lower():
               with open(RES_DIR+"pos_sequence.fasta", "a") as fasta_file:
                    fasta_file.write(capitalise(seq))

# create positive dataset with annotations
with open(RES_DIR + "pos_sequence.fasta", "r") as file:
    data = file.read()
    seq_list = data.split(">")

    for i in tqdm(range(len(seq_list))):
        pos_seq_list = seq_list[i].strip().split("\n")
        if len(pos_seq_list) == 2:
            pdb_id = pos_seq_list[0][0:6].strip()
            #position_list = [int(x[0].split('+')[0]) for x in pos_dict[pdb_id]]  # in files positions start from 1
            position_list = [int(x[0]) for x in pos_dict[pdb_id]]  # in files positions start from 1
            position_list.sort()
            seq = insert(pos_seq_list[1].strip(),position_list,"#")
            print(pdb_id , position_list, seq)
            with open(RES_DIR + "pos_annotated_sequence.fasta", "a") as fasta_file:
                fasta_file.write(">"+pdb_id+"\n"+seq)
            #print(pos_seq_list[1])
        # pdb_id = pos_seq_list[0][0:6]
        # seq = pos_seq_list[1]
        # print(seq)

        # print(len(pos_seq_list))