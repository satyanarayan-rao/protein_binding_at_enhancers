import os
import sys
import re
from collections import defaultdict
import pickle
cluster_fp = open(sys.argv[1])
footprint_length_dict = pickle.load(open(sys.argv[2], "rb"))
out_fp = open(sys.argv[3], "w")
out_fp_with_per_orange = open(sys.argv[4], "w")
cl_id_len_dict = defaultdict(list)
for line in cluster_fp: 
    annotation = line [0:line.find("\t")]
    footprint_str = line[line.find("\t") + 1:len(line) - 1 ]
    per_orange = round(footprint_str.count("F")*100/(len(footprint_str)), 3) 
    cl_id = annotation[annotation.find("#") + 1: len (annotation)]
    chr_loc = annotation[0: annotation.find("#")]
    if len(footprint_length_dict[chr_loc]) == 0: 
        cl_id_len_dict[cl_id].extend([0]) 
        to_write = annotation + "\t" + str(0) + "\t" + str(0)
        out_fp_with_per_orange.write(to_write + "\n")
    else:
        cl_id_len_dict[cl_id].extend(footprint_length_dict[chr_loc]) 
        for l in footprint_length_dict[chr_loc]: 
            to_write = annotation + "\t" + str(per_orange) + "\t" + str(l)
            out_fp_with_per_orange.write(to_write + "\n")

all_keys = list(cl_id_len_dict.keys()) 
all_keys.sort() 
cl = 1 
out_fp.write("cl_id\tfootprint_len\n")
for k in all_keys:
     for j in cl_id_len_dict[k]:
         to_write = str(cl) + "\t" + str(j)
         out_fp.write(to_write + "\n")
     cl += 1 
out_fp.close() 
cluster_fp.close()
out_fp_with_per_orange.close()
