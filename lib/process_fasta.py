import sys
import pandas as pd
import numpy as np
import os
from sklearn.model_selection import KFold,RepeatedKFold

RANDOM = 8


prefix = ""
fasta_dir = "../../data/resource/K-Fold/"
fasta_filename = "../../data/resource/annotated_sequence.fasta"


def write_dict_fastafile(filename = "", fasta_dict = {}):
    """write_dict 2fastafile
    Args:
        filename (TYPE, optional):
        fasta_dict (dict, optional):
    """
    file_write_obj = open(fasta_dir+"annotated_sequence.fasta", 'w')
    for k,v in fasta_dict.items():
        file_write_obj.writelines(k)
        file_write_obj.write('\n')
        file_write_obj.writelines(v)
        file_write_obj.write('\n')
    file_write_obj.close()

def read_fasta(fasta_file):
    """read fasta to dict
    convert X to A
    Args:
        fasta_file (TYPE):fasta  with annotation #
    Returns:
        fasta_dict {}, annotated
        idlist: []
    """
    fp = open(fasta_file)
    lines = fp.readlines()

    fasta_dict = {} #record seq for one id
    positive_dict={} #record positive positions for one id
    idlist=[] #record id list sorted
    gene_id = ""
    for line in lines:
        if line[0] == '>':
            if gene_id != "":
                fasta_dict[gene_id] = seq
                idlist.append(gene_id)
            seq = ""
            gene_id = line.strip('\n') #  line.split('|')[1] all in > need to be id
        else:
            seq += line.strip('\n')
            seq = seq.replace('X','A')
            # seq = seq.replace('#','')
    fasta_dict[gene_id] = seq #last seq need to be record
    idlist.append(gene_id)

    return fasta_dict, idlist

def replace_anno(fasta_dict):
    """replace_annotaion # with ""
    counting number of # of each id
    Args:
        fasta_dict (TYPE): Description
    Returns:
        fasta_dict_,
        id_len_dict,
        id_annoNum_dict
    """
    fasta_dict_ = {}
    id_len_dict = {}
    id_annoNum_dict = {}
    for k,v in fasta_dict.items():
        id_annoNum_dict[k] = v.count("#")
        fasta_dict_[k] = v.replace('#','')
        id_len_dict[k] = len(v)
    print("id_len_dict, id_annoNum_dict")
    print(id_len_dict.values(), id_annoNum_dict.values())
    return fasta_dict_, id_len_dict, id_annoNum_dict


def cdhit(th = 0.4, word_size_n = 2, input_filename = "", output_filename = ""):
    """
    Choose of word size:
    -n 5 for thresholds 0.7 ~ 1.0
    -n 4 for thresholds 0.6 ~ 0.7
    -n 3 for thresholds 0.5 ~ 0.6
    -n 2 for thresholds 0.4 ~ 0.5 # default
    https://github.com/weizhongli/cdhit/wiki/3.-User%27s-Guide#CDHIT
    Args:
        th (float, optional):
        word_size_n (int, optional):
        input_filename (str, optional):
        output_filename (str, optional):
    """
    if output_filename == "": output_filename = input_filename+str(th)
    os.system("cd-hit -i "+input_filename+" -o "+output_filename+" -c "+str(th)+" -n "+str(word_size_n)+" -M 16000 -T 8 -d 0")

def read_cdhit_clstrfile(clstrfilename = fasta_dir+"annotated_sequence.fasta0.4.clstr"):
    """Summary
    Args:
        clstrfilename (str, optional): ">Cluster 3"
    Returns:
        cls_dict
    """
    cls_dict = {}
    fp = open(clstrfilename)
    lines = fp.readlines()

    positive_dict={} #record positive positions for one id
    cls_id = ""
    seq = []
    for line in lines:
        if line.startswith('>'):
            if cls_id != "":
                cls_dict[cls_id] = seq
            seq = []
            cls_id = int(line.strip('\n').split()[1]) #  line.split('|')[1] all in > need to be id
        else:
            fasta_id = line.split('\t')[1].split()[1]
            fasta_id = fasta_id[:-3]
            print(fasta_id)
            seq.append(fasta_id)

    cls_dict[cls_id] = seq #last seq need to be record
    print(cls_dict)
    return cls_dict

def vali_test_cls_pick(test_index, cls_dict,  fasta_dict_anno, id_len_dict = None, id_annoNum_dict = None, pick_method = "annoNum"):
    """vali and test extraction for each fold, pick 1 sequence
    test_index: [int...], inedx in cluster
    pick_method = "annoNum" or "len" or None
        annoNum: pick the sequence with most #
        len: pick the longest sequence
        None: pick all sequences
    Return:
    vali_dict[fastaid]={seq}
    """
    vali_dict = {}

    for clsid in test_index:
        fastaid_arr = np.array([fastaid for fastaid in cls_dict[clsid]])
        annoNum_arr = np.array([id_annoNum_dict[fastaid] for fastaid in cls_dict[clsid]])
        len_arr = np.array([id_len_dict[fastaid] for fastaid in cls_dict[clsid]])

        if pick_method != None:
            max_id_in_arr = np.argmax(eval(pick_method+"_arr"))
            fastaid_pick = fastaid_arr[max_id_in_arr]
            vali_dict[fastaid_pick] = fasta_dict_anno[fastaid_pick]
        else:
            for i in fastaid_arr:
                vali_dict.update(fasta_dict_anno[i])


    return vali_dict

def trn_tst_KFold(cls_dict, fasta_dict_anno, id_len_dict = None, id_annoNum_dict = None, fold_num = 9, n_repeats=1):
    """write fold files by train test and validation
    Args:
        X []: cls_ids
        y []: cls_fastaid
    """

    cls_ids = np.array(cls_dict.keys())
    cls_fastaid = np.array(cls_dict.values())

    kf = RepeatedKFold(n_splits=fold_num, n_repeats=n_repeats, random_state=RANDOM)
    i = 0
    for train_index, test_index in kf.split(cls_ids):
        print("TRAIN:", train_index, "TEST:", test_index)
        train_dict_fold = {}
        for clsid in train_index:
            for fastaid in cls_dict[clsid]:
                train_dict_fold[fastaid] = fasta_dict_anno[fastaid]

        write_dict_fastafile(filename = fasta_filename+"_training_annotated_"+str(i)+".fasta", fasta_dict = train_dict_fold)

        test_dict = vali_test_cls_pick(test_index[: len(test_index)/2], cls_dict, fasta_dict_anno, id_len_dict = id_len_dict, id_annoNum_dict = id_annoNum_dict,pick_method=None)
        vali_dict = vali_test_cls_pick(test_index[len(test_index)/2 :], cls_dict, fasta_dict_anno, id_len_dict = id_len_dict, id_annoNum_dict = id_annoNum_dict,pick_method=None)
        write_dict_fastafile(filename = fasta_filename+"_testing_annotated_"+str(i)+".fasta", fasta_dict = test_dict)
        write_dict_fastafile(filename = fasta_filename+"_validation_annotated_"+str(i)+".fasta", fasta_dict = vali_dict)

        i += 1
    return train_dict_fold, vali_dict, test_dict



def processing_anno_fasta():

    fasta_dict_anno, _ = read_fasta(fasta_filename)
    fasta_dict_new, id_len_dict, id_annoNum_dict = replace_anno(fasta_dict_anno)
    write_dict_fastafile(filename = fasta_filename, fasta_dict = fasta_dict_new)
    cdhit(input_filename = prefix+fasta_filename)

    cls_dict = read_cdhit_clstrfile(clstrfilename = prefix+fasta_filename+"0.4.clstr")

    trn_tst_KFold(cls_dict, fasta_dict_anno, id_len_dict, id_annoNum_dict)

if __name__ == '__main__':
    processing_anno_fasta()