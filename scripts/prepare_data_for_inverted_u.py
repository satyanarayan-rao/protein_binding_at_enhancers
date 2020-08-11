import os
import sys
import re
import pickle
from collections import defaultdict

cluster_fp  = open(sys.argv[1])
footprint_length_per_bp_dict = pickle.load(open(sys.argv[2], "rb")) 
out_fp = open(sys.argv[3], "w")
#print(list(footprint_length_per_bp_dict.keys)[0:5])

for line in cluster_fp:
    annotation = line [0:line.find("\t")]
    chr_loc = annotation[0: annotation.find("#")] 
    footprint_len_per_bp = "\t".join(map(str, footprint_length_per_bp_dict[chr_loc])) 
    to_write = annotation + "\t" + footprint_len_per_bp
    out_fp.write(to_write + "\n")
out_fp.close() 
cluster_fp.close() 
