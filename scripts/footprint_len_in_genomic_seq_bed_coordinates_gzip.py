import os 
import sys
import re
import gzip
inp_fp = open(sys.argv[1])
out_fp = gzip.open(sys.argv[2], "wb")
from itertools import combinations
from collections import defaultdict 
# find CpG and GpC in genomic sequence and calculate distance amongst them
# example sequence
#             0123456789012345678901234567890 - total 31 bp   (-15 to 15)
# peak_2237`1	TAGTCAAGCGAAGTTCGGCGACAACAGCAGT
# potential C --------|------|--|--------|---
# Footprint 1 ---------FFFFFFFFF-------------
# Footprint 2 ---------FFFFFFFFFFFFFFFFFF---- # result by either of 2 middle unmethylated cytosines
# Fooptrint 3 ----------------FFFFFFFFFFF----
#                     
# there are only 3 Cs found in either CpG or GpC context in this case (which is required Methylated - Unmethylated - Methylated)


for line in inp_fp:
    line_items = line.strip().split()
    peak_id = "\t".join(line_items[0:3])
    read_id = line_items[-1]
    start = int(line_items[1])
    sequence = line_items[3]
    dinuc_loc = []
    for idx in range (len(sequence) - 1):
        if sequence[idx:idx+2] == "CG":
            dinuc_loc.append(idx) # locating | for CG shown in example
        elif sequence[idx:idx+2] == "GC": 
            dinuc_loc.append(idx + 1) # locating | for GC shown in example
            
    # since the location are alreay sorted in nature - taking combination with r = 3 will ensure footprint call and length will be thid minus first minus 1
    if len(dinuc_loc)>2:
        all_comb = list(combinations(dinuc_loc, 3))
        footprint_start_dict = defaultdict(list)
        for j in all_comb:
            dist = j[2] - j[0] - 1  # excluding both pipes hence -1  
            footprint_start_dict[j[0]].append(dist)
        for k in footprint_start_dict:
            footprint_list = set(footprint_start_dict[k])
            for l in footprint_list:
                to_write = line_items[0] + "\t" +  str(start + k) + "\t" + str (start + k + l) + "\t" +  str(l) + "\t" + read_id
                out_fp.write(bytes(to_write+"\n", encoding = "ascii"))
    else:
        to_write = peak_id + "\t" + "NA" + "\t" + read_id
        out_fp.write(bytes(to_write + "\n", encoding = "ascii"))
out_fp.close()
inp_fp.close()
