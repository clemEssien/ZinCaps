import os
import shutil
from glob import glob
from collections import defaultdict
import re
from tqdm import tqdm
import textwrap
import requests
from lxml import html
from lib import process_fasta
from lib import library
DATASET_DIR = "../../data/dataset/"
RES_DIR = "../../data/resource/"
DSSP_DIR_OLD = "../../data/resource/old_dssp"
DSSP_DIR_NEW = "../../data/resource/dssp"
DSSP_DIR_IND =  "../../data/resource/dssp_independent_dataset"
positive_pids = []
negative_pids = []
temp_pos_list = []
temp_neg_list = []
pos_dict = defaultdict(list)
neg_dict = defaultdict(list)

def move_file_to_folder(pdb_id, folder_old, folder_new):
    for filename in os.listdir(folder_old):
        filename = filename.split('.dssp')[0]
        if filename.strip().lower() == pdb_id.strip().lower():
            shutil.move(folder_old+"/"+filename+".dssp", folder_new+"/"+filename+".dssp")

def write_pos_annotations(pdb_id,seq,filename):
    with open(RES_DIR + filename, "a") as fasta_file:
        fasta_file.write(">"+pdb_id+"\n"+textwrap.fill(seq,60)+"\n")


def annotate(file_id,pos_dict):
    with open(DSSP_DIR_NEW + "/" + file_id[0:4].lower() + ".dssp", "r") as file:
        line_no = 0
        seq = ""
        for line in file:
            line_no += 1
            if line_no > 28 and line.split()[len(line.split()) - 1].strip() == file_id[5:6].strip().upper():
                seq += line.split()[3]
                if file_id in positive_pids:
                    pos = [(x[0]) for x in pos_dict[file_id]]
                    if line.split()[1] in pos:
                        seq += "#"

    return seq

def create_dssp(acc_id):
    API_ENDPOINT ='https://mrs.cmbi.umcn.nl/search?db=dssp&q='+acc_id

    # sending post request and saving response as response object
    response = requests.get(url = API_ENDPOINT, stream=True)
    tree = html.fromstring(response.content)

    dssp = tree.xpath('//pre[@id="entrytext"]/text()')

    with open(DSSP_DIR_NEW+'/'+acc_id[0:4].lower()+'.dssp', 'w') as f:
        for ds in dssp:
            f.write(ds)
    print("successfully created file : "+acc_id[0:4].lower()+".dssp")

def parse_dssp_file(file_id, pos_dict,filename):
    seq = annotate(file_id, pos_dict)
    write_pos_annotations(file_id, seq, filename)

# loop through both positive and negative dataset files
# to create a dictionary comprising of accesion-id as key and the position(s) as values
# positions follow annotations in the DSSP
for files in glob(DATASET_DIR+"*"):
    filename = os.path.split(files)[1]
    file = open(DATASET_DIR+filename,"r")
    for f in file:
        if "positive" in filename :
            id = f[0:6] # retrieve id
            positive_pids.append(id);
            if positive_pids[len(positive_pids)-1] == id :
                temp_pos_list.append(f[11:14].strip()) # get the positions
            temp_pos_list.sort()
            pos_dict[id].append(temp_pos_list) # creates dictionary have id as keys and list of annotations as value
            temp_pos_list = []
        elif "negative" in filename:
            id = f[0:6]
            negative_pids.append(id);
            if negative_pids[len(negative_pids)-1] == id:
                temp_neg_list.append(f[11:14].strip())
            temp_neg_list.sort()
            neg_dict[id].append(temp_neg_list)
            temp_neg_list = []

positive_pids = list(set(positive_pids))
negative_pids = list(set(negative_pids))
combined_pids = list(set(positive_pids +negative_pids))
print("total ids: " ,len(combined_pids))

# use the ids from the dictionary of combined annotations to open the corresponding dssp file to generate the sequence
for i in tqdm(range(len(combined_pids))):
    parse_dssp_file(combined_pids[i], pos_dict,"annotated_sequence.fasta")

#fasta_dict, idlist = library.read_fasta(RES_DIR+"sequence.fasta")
#print(len(fasta_dict))

# selecting only sequences that have at least one annotated site. But I'm avoiding this now to improve the
# specificity
# fasta_dict_new = {}
# for k,v in fasta_dict.items():
#     #print(k,v)
#     if "#" in fasta_dict[k]:
#         fasta_dict_new[k] = fasta_dict[k]



# library.write_dict_fastafile(RES_DIR+"annotated_sequence", fasta_dict_new)
#
# process_fasta.processing_anno_fasta()