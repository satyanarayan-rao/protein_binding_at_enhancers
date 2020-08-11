import os
import sys
import re
from collections import defaultdict
cluster_fp = open(sys.argv[1])
footprint_length_fp = open(sys.argv[2])
lflank = int(sys.argv[3])
rflank = int(sys.argv[4])
out_fp = open (sys.argv[5], "w")


# first create a dict of all length in without matrix form
flen_list_dict = defaultdict(lambda : defaultdict(list)) 
for line in footprint_length_fp:
    d_loc = [m.start() for m in re.finditer("\t", line)]
    k = "`".join ( [ line[d_loc[2] + 1: d_loc[3]],
                     line[d_loc[8] + 1: d_loc[9]] ])
    #print (k)
    #sys.exit(-1)
    footprint_start_on_read = int (line[d_loc[-2] + 1: d_loc[-1]].split("-")[0]) + 1 
    footprint_end_on_read = int (line[d_loc[-2] + 1: d_loc[-1]].split("-")[1]) + 1  
    footprint_length = int(line[d_loc[-1] +1: len(line)]) 
    relative_pos = int(line[d_loc[0] + 1: d_loc[1]]) - int(line[d_loc[6] + 1: d_loc[7]])  -  footprint_start_on_read
    flen_list_dict[k][relative_pos].append(footprint_length)

flist_dict = defaultdict(list)
for line in cluster_fp: 
    line_items = line.strip().split()
    t = line_items[0].split ("`") 
    k = t[0] + "`" + t[1]
    cl_id = int(line_items[0].split("#")[1])
    for kk in flen_list_dict[k].keys(): 
        if (kk >= lflank) and (kk<=rflank):
            flist_dict[cl_id].extend(flen_list_dict[k][kk])
        else: 
            for j in flen_list_dict[k][kk]: 
                if j + kk >= lflank:
                    flist_dict[cl_id].append(j)

all_keys = list(flist_dict.keys()) 
all_keys.sort() 

cl_id = 1  
out_fp.write("cl_id\tfootprint_len\n")
for k in all_keys:
    for j in flist_dict[k]:
        to_write = "\t".join([str(cl_id), str(j)])
        out_fp.write(to_write + "\n")
    cl_id += 1 

out_fp.close()        
