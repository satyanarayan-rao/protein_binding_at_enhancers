import os
import sys
import pickle
import re

inp_fp = open(sys.argv[1])
mapped_read_dict = pickle.load(open(sys.argv[2], "rb"))
peak_id = sys.argv[3]
s_peak_loc = int(sys.argv[4])
strand = sys.argv[5]
sam_flag = sys.argv[6]
out_fp = open(sys.argv[7], "w")
lf = int(sys.argv[8])
rf = int(sys.argv[9])
for line in inp_fp:
    # example line
    # chr3L:632728-632728^10`SRR3133326.1889480_1889480/1_overlapping`99~147	MMMMMMMMMMMMMMMMMMMM.........................................................................................................................................................................................................................................................................................
    read_id = line[line.find("`") + 1: line.find("\t")] 
    read_start = mapped_read_dict[read_id]["start"]
    read_end = mapped_read_dict[read_id]["end"]
    fp_center =  int(line[line.find(":") + 1: line.find("-")])
    
    if (strand == "+") or (strand == "."):

        if (fp_center - lf >= read_start) and (fp_center + s_peak_loc + rf <= read_end):
            out_fp.write(line)
        else: 
            continue 
    elif strand == "-": 
        if (fp_center + lf <= read_end) and (fp_center - s_peak_loc - rf >= read_start):
            out_fp.write(line)
        else: 
            continue 
        
out_fp.close()
