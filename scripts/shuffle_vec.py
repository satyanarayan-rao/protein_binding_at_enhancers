import os
import sys
import re
from collections import defaultdict, OrderedDict
import potential_cytosines
import numpy as np
import random

def find_idx_all_potential_cytosines(genomic_seq, sam_flag, cytosine_dict):
    """
    genomic_seq = AGAAGAGCTGTATGTGCGTGCATCTCCTCGACCAGTCAGTCTTCTCTGCTTTTCTTCTAGATTCTAGACGGGGGGCCTCTCCTCGTTCTCTCCCTCTCACCATTCTGGTTTGCTTTTGGCTGCCAGCAGGGAAAAAAAGTTGTGTCTTTCTTGCCATGCCGCCTCCTCTTTGTTCCCTTCTTACCCACTGCTTTTCGTGCACTTTTGTACACGGTTTTGCTTTTTTTATTTTGTTATTCTCGTCTCGCCTCTCTCTCGCTCCTGCCTCTCTCCCTCTCTT 
                  .h..h.H..x...h.H.Z.H.........Z....x...x........X...........h......h..zxhhhh.........z.....................xh...h.....hh..x...x..xhh.......h..h.h........h....h..z..........h.................X......Z.H.......h.....Zx....H.............h........Z....Z..........z.....X................
    sam_flag = "83~163"
    cytosine_dict ~ see scripts/potential_cytosines.py
    """
    idx_list = [] 
    seq_len = len(genomic_seq)
    for i in range(seq_len - 3):
        if cytosine_dict[sam_flag][genomic_seq[i:i+3]]: 
            first = genomic_seq[i:i+2]
            second = genomic_seq[i+1:i+3]
            if sam_flag == "99~147":
                idx_list.append(i)
                if first == "CG":
                    idx_list.append(i)
                elif second == "CG":
                    idx_list.append(i+1)
            else:
                idx_list.append(i+2) 
                if first == "CG":
                    idx_list.append(i+1)
                elif second == "CG":
                    idx_list.append(i+2)
                
            
        elif cytosine_dict[sam_flag][genomic_seq[i:i+2]]:
            if sam_flag == "99~147":
                idx_list.append(i)
            else:
                idx_list.append(i+1)            
        else:
            continue 
    if cytosine_dict[sam_flag][genomic_seq[seq_len - 2:seq_len]]:
        idx_list.append(seq_len - 2) 
    unique_idx_list = list(np.unique(idx_list)) 
    return unique_idx_list 
def get_all_non_dot_loc(mvec):
    result = []
    methylated_count = 0
    unmethylated_count = 0 
    for i in range(len(mvec)):
        if mvec[i]!=".":
            result.append(i)
            if mvec[i] == mvec[i].upper(): 
                methylated_count +=1
            else: 
                unmethylated_count +=1 
    return result, methylated_count, unmethylated_count

genomic_fp = open(sys.argv[1])
mvec_fp = open(sys.argv[2])
shuf_fp = open(sys.argv[3], "w")
verbose_fp = open (sys.argv[4], "w")
random.seed(int(sys.argv[5]))

cdict = potential_cytosines.get_potential_cytosine_dict()

for g_line, m_line in zip(genomic_fp,mvec_fp): 
    #chr2L	11806729	11806730	chr2L:11806729-11806729^20	.	+	chr2L	11806644	11806923	SRR3133326.203452_203452/1_overlapping`83~163	.	AGAAGAGCTGTATGTGCGTGCATCTCCTCGACCAGTCAGTCTTCTCTGCTTTTCTTCTAGATTCTAGACGGGGGGCCTCTCCTCGTTCTCTCCCTCTCACCATTCTGGTTTGCTTTTGGCTGCCAGCAGGGAAAAAAAGTTGTGTCTTTCTTGCCATGCCGCCTCCTCTTTGTTCCCTTCTTACCCACTGCTTTTCGTGCACTTTTGTACACGGTTTTGCTTTTTTTATTTTGTTATTCTCGTCTCGCCTCTCTCTCGCTCCTGCCTCTCTCCCTCTCTT
    d_loc = [m.start() for m in re.finditer("\t", g_line)]
    sam_flag = g_line[d_loc[8] + 1: d_loc[9]].split("`")[-1]
    genomic_seq = g_line[d_loc[-1] + 1: len(g_line) - 1] # -1 is to avoid `\n`
    print ([sam_flag, genomic_seq])
    all_potential_cytosines =  find_idx_all_potential_cytosines(genomic_seq, sam_flag, cdict) 
    #print(all_potential_cytosines)
    #print (m_line)
    # record all marked positions 
    mvec = m_line[d_loc[-1] + 1: len(g_line) - 1]
    res, methylated_count, unmethylated_count = get_all_non_dot_loc(mvec)
    # First find common between mvec and potential
    genomic_set = set(all_potential_cytosines) 
    mvec_set = set(res)
    # common 
    set_dict = OrderedDict()
    set_dict["common"] = sorted(list(set.intersection(genomic_set, mvec_set)))
    set_dict["unique_in_genomic"] = sorted(list(set.difference(genomic_set, mvec_set)))
    set_dict["unique_in_mvec"] = sorted(list(set.difference(mvec_set, genomic_set)))
    to_write = ""
    to_append = ""
    for k in set_dict:
        if len(set_dict[k]) == 0:
            to_append = k + "@" + "NA"
        else:
            to_append = k + "@" + ",".join(map(str, set_dict[k]))
        if to_write == "":
            to_write = to_append
        else:
            to_write = "|".join([to_write, to_append])
    to_write = g_line[0:d_loc[10]] + "\t" + to_write 
    verbose_fp.write(to_write + "\n")
    # shuffle the positions
    mvec_list = list(mvec)
    res_len = len(res)
    shuffle_methylated = random.sample(range(res_len), methylated_count)
    for i in shuffle_methylated: 
        mvec_list[res[i]] = "Z"
    for i in range(res_len):
        if i not in shuffle_methylated:
            mvec_list[res[i]] = "z" 
    shuf_vec = "".join(mvec_list)
    shuf_write = g_line[0:d_loc[10]] + "\t" + shuf_vec
    shuf_fp.write(shuf_write + "\n") 
        
verbose_fp.close()
shuf_fp.close()
       
    
