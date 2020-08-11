import sys
import pickle
import re 
import numpy as np

inp_fp = open(sys.argv[1])
peak_dict = pickle.load(open(sys.argv[2], "rb"))
out_fp = open(sys.argv[3], "w")

for line in inp_fp: 
    # example line: head -1 all_peak_pos
    # peak_5354	2R	POS:15519001	DIST:0	POS:15518715	DIST:286	POS:15518868	DIST:133
    all_pos = re.findall(r"((?<=POS:)(\d+))", line)
    all_pos_l =  [int (i[1]) for i  in all_pos]
    all_dist = re.findall(r"((?<=DIST:)(\d+))", line)
    all_dist_l =  [int (i[1]) for i  in all_dist]
    first_tab = line.find("\t")
    second_tab = line[first_tab + 1:].find ("\t")
    peak_id = line[0:first_tab]
    peak_chr = "chr" + line[first_tab+1: first_tab + second_tab + 1] 
    
    sort_idx = np.argsort(all_pos_l)
    sorted_all_pos_l = [all_pos_l[i] for i in sort_idx]
    sorted_all_dist_l = [sorted_all_pos_l[i] - sorted_all_pos_l[0] for i in range(len(sorted_all_pos_l))]
    cl_id = None
    if peak_id in peak_dict: 
    	  cl_id = peak_dict[peak_id]["cl_id"]
    else: 
        cl_id = "100"
    for pos, dist in zip(sorted_all_pos_l, sorted_all_dist_l):
        to_write = "\t".join([peak_chr, str(pos), str(pos), str(dist), peak_id, cl_id, "+" ]) 
        out_fp.write(to_write + "\n")
        
