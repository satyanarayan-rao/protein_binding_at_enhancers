import os
import sys
import pickle
import re

input_fp = open(sys.argv[1])
out_dict = open(sys.argv[2], "wb")
read_methylation_dict = {}
cnt = 1 
for line in input_fp:
    #print (cnt)
    line = line.strip() 
    len_line = len(line)
    tmp_dict = {}
    all_tabs = [m.start() for m in re.finditer("\t", line)] 
    #print (all_tabs)
    chrom_loc = int(line[all_tabs[1]: all_tabs[2]])
    methylation_vec = line[all_tabs[2] + 1: len_line]
    for i in range(len(methylation_vec)): 
        m_status = methylation_vec[i]
        if m_status == ".":
            continue
        else: 
            
            dict_key = ":".join([
                        line[0:all_tabs[0]], # read name 
                        line[all_tabs[0]+1:all_tabs[1]],
                        str(chrom_loc + i)])
            read_methylation_dict[dict_key] = m_status
    cnt = cnt + 1 
pickle.dump(read_methylation_dict, out_dict) 
out_dict.close() 
input_fp.close()
