import os
import sys
import re
from collections import defaultdict 

inp_fp = open(sys.argv[1])
out_fp = open (sys.argv[2], "w")
label_dict = {
    "0" : "occup_D_D",
    "1" : "occup_D_T",
    "2" : "occup_D_N",
    "3" : "occup_T_D",
    "4" : "occup_T_T",
    "5" : "occup_T_N",
    "6" : "occup_N_D",
    "7" : "occup_N_T",
    "8" : "occup_N_N"
}


#header = "file_name\ttotal_reads\toccup_naked\toccup_tf\toccup_nuc"
header = "\t".join(["file_name", "total_reads", "\t".join([label_dict[str(i)] for i in range(len(label_dict))])]) 

out_fp.write (header + "\n") 
for line in inp_fp: # each line is a file name
    fname = line.strip()
    fp = open(fname) 
    total_reads = 0
    occup_count_dict = defaultdict(lambda : 0)
    for f_line in fp:
        total_reads +=1 
        binding_id = f_line[f_line.find("#")+1: f_line.find("\t")] 
        occup_count_dict[label_dict[binding_id]] +=1 
    fp.close()     
    percent_vec = []
    for k in ["occup_D_D", "occup_D_T", "occup_D_N","occup_T_D", "occup_T_T", "occup_T_N","occup_N_D", "occup_N_T", "occup_N_N"]: 
        percent_vec.append (round((occup_count_dict[k]/total_reads)*100,2))
    to_write = "\t".join(map(str,percent_vec))
    to_write = fname + "\t" + str(total_reads) +"\t" + to_write
    out_fp.write(to_write + "\n") 

out_fp.close()
inp_fp.close()

