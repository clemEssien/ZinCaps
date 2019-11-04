import re
import os
import textwrap
def isEmpty(list):
    if not list:
        return True
    else:
        return False

def print_list(list):
	for i in range(len(list)):
		print(list[i]);

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
    sequence.append(array[i]+1) #seems the regphos position annotation counts from zero
    new_string += string[sequence[i]:sequence[i+1]]+xter
  new_string += string[sequence[-1]:]+"\n"
  return new_string

def capitalise(str):
    return ">" + str[0:6].upper() + str[6:len(str)]

def fileList(folder_old):
    files = []
    for filename in os.listdir(folder_old):
        filename = filename.split('.dssp')[0]
        files.append(filename)
    return files

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
            gene_id = line.strip('\n').strip('\r') #  line.split('|')[1] all in > need to be id
        else:
            seq += line.strip('\n')
            #seq = seq.replace('X','A')
            # seq = seq.replace('#','')
    fasta_dict[gene_id] = seq #last seq need to be record
    idlist.append(gene_id)
    return fasta_dict, idlist

def write_dict_fastafile(filename,fasta_dict):
    """write_dict 2fastafile
    Args:
        filename (TYPE, optional):
        fasta_dict (dict, optional):
    """
    file_write_obj = open(filename, 'w')
    for k,v in fasta_dict.items():
        file_write_obj.writelines(k)
        file_write_obj.write('\n')
        file_write_obj.writelines(textwrap.fill(v,60))
        file_write_obj.write('\n')
    file_write_obj.close()