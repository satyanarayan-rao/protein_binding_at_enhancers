import os
import sys
import re
from collections import defaultdict 

inp_fp = open(sys.argv[1])
out_fp = open (sys.argv[2], "w")
label_dict = {
    "0" : "Naked-DNA",
    "1" : "TF",
    "2" : "Nuc"
    "3" : "MostlyNuc"
}


header = "file_name\ttotal_reads\toccup_naked\toccup_tf\toccup_nuc"
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
    for k in ["Naked-DNA", "TF", "Nuc", "MostlyNuc"]: 
        percent_vec.append (round((occup_count_dict[k]/total_reads)*100,2))
    to_write = "\t".join(map(str,percent_vec))
    to_write = fname + "\t" + str(total_reads) +"\t" + to_write
    out_fp.write(to_write + "\n") 

out_fp.close()
inp_fp.close()

