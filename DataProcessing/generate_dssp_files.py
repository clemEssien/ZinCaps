import requests
import os
import textwrap
from lxml import html
from glob import glob
from collections import defaultdict
from tqdm import tqdm

RES_DIR = "../../data/resource/"
DATASET_DIR = "../../data/dataset/"
DSSP_DIR_NEW = "../../data/resource/dssp/"
pos_dict = defaultdict(list)
temp_pos_list = []

def write_pos_annotations(pdb_id,seq,filename):
    with open(DSSP_DIR_NEW + filename, "a") as fasta_file:
        fasta_file.write(">"+pdb_id+"\n"+textwrap.fill(seq,60)+"\n")


def parse_dssp_file(file_id, pos_dict,filename,acc_ids):
    with open(DSSP_DIR_NEW + file_id[0:4].lower() + ".dssp", "r") as file:
        line_no = 0
        seq = ""
        for line in file:
            line_no += 1
            if line_no > 28 and line.split()[len(line.split())-1].strip() == file_id[5:6].strip().upper():
                seq += line.split()[3]
                if file_id in acc_ids:
                    pos = [(x[0]) for x in pos_dict[file_id]]
                    if line.split()[1] in pos:
                        seq += "#"
    write_pos_annotations(file_id, seq, filename)


def create_dssp(acc_id):
    API_ENDPOINT ='https://mrs.cmbi.umcn.nl/search?db=dssp&q='+acc_id

    # sending post request and saving response as response object
    response = requests.get(url = API_ENDPOINT, stream=True)
    tree = html.fromstring(response.content)

    dssp = tree.xpath('//pre[@id="entrytext"]/text()')


    with open(DSSP_DIR_NEW+acc_id.lower()+'.dssp', 'w') as f:
        for ds in dssp:
            f.write(ds)
    print("successfully created file : "+acc_id+".dssp")


acc_ids = []
for files in glob(DATASET_DIR+"*"):
    filename = os.path.split(files)[1]
    file = open(DATASET_DIR+filename,"r")
    for f in file:
        if "zhao" in filename :
            id = f[0:6]
            create_dssp(id[0:4])
            acc_ids.append(id)
            if acc_ids[len(acc_ids) - 1] == id:
                temp_pos_list.append(f[11:14].strip())  # get the positions
                #print(id,f[11:14])
            temp_pos_list.sort()
            pos_dict[id].append(temp_pos_list)  # creates dictionary have id as keys and list of annotations as value
            temp_pos_list = []
            # if(id not in acc_ids):
            #     acc_ids.append(id)

for i in tqdm(range(len(acc_ids))):
    parse_dssp_file(acc_ids[i], pos_dict,"ind_sequence.fasta",acc_ids)

#print(acc_ids)
#print(len(acc_ids))
# for acc_id in acc_ids:
#     print("creating "+acc_id[0:4]+ " file.....")
#     create_dssp(acc_id[0:4])