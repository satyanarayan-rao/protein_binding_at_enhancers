import sys
from collections import defaultdict
import re

footprint_fp = open(sys.argv[1])
mvec_fp = open(sys.argv[2])
out_fp = open(sys.argv[3], "w")

read_id_label_dict = {} 
for line in footprint_fp: 
    read_id_and_label = line[0:line.find("\t")].split("#")
    read_id_label_dict[read_id_and_label[0]] = read_id_and_label[1]

for line in mvec_fp:
    d_loc = line.find("\t")
    read_id = line[0:d_loc] 
    if read_id in read_id_label_dict: 
        out_fp.write(read_id + "#" + read_id_label_dict[read_id] + "\t" + line[d_loc + 1:]) 
    else: 
        continue

out_fp.close()
