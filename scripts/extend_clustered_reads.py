import os
import sys
import pickle
import re

input_fp = open(sys.argv[1])
footprint_pkl = pickle.load(open (sys.argv[2], "rb"))
out_fp = open(sys.argv[3], "w")
original_l = int(sys.argv[4])
original_r = int(sys.argv[5])
extend_l = int(sys.argv[6])
extend_r = int(sys.argv[7])
for line in input_fp: 
    # get the first element of line before "#"
    key_with_hash = line[0:line.find("\t")]
    footprint_str = line[line.find("\t") + 1: len(line) - 1]
    key = key_with_hash.split("#")[0]
    if key in footprint_pkl: 
        extended_footprint = footprint_pkl[key]
        out_fp.write(key_with_hash + "\t" + extended_footprint + "\n")
    else: 
        missing_filled = "M"*(extend_l - original_l) + footprint_str + "M"*(extend_r - original_r) 
        out_fp.write(key_with_hash + "\t" + missing_filled + "\n")

out_fp.close() 
input_fp.close() 
