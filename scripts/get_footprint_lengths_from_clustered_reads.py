import os
import sys
from collections import defaultdict
# sys.argv[1]: clusted file
# chr3R:19134673-19134673^12`SRR3133326.927598_927598/1_overlapping`99~147#63     ..............................................................FFF.....FFFFFFFFFFFFF....................................................................

flen_dict_list = defaultdict(list)
in_fp = open(sys.argv[1]) 
out_fp = open (sys.argv[2], "w")


for line in in_fp:
    line_items = line.strip().split()
    cl_id = int(line_items[0].split ("#")[-1])
    flen = 0
    for c in line_items[1]:
        if c == "." and flen == 0:
            continue
        elif c == "." and flen !=0: 
            flen_dict_list[cl_id].append(flen)
            flen = 0
        else:
            flen +=1
all_keys = list(flen_dict_list.keys()) 
all_keys.sort() 
print(all_keys)

cl_id = 1 
for k in all_keys:
    for j in flen_dict_list[k]:
        to_write = "\t".join([str(cl_id), str(j)])
        out_fp.write(to_write + "\n")
    cl_id += 1 

out_fp.close()
