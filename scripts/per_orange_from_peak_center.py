import os
import sys
import re
from collections import defaultdict
import pickle
# Look at the rule get_real_footprint_lengths_and_oraange_on_read

inp_fp = open(sys.argv[1])
footprint_length_dict = pickle.load(open(sys.argv[2], "rb"))
out_fp_per_orange_and_flen = open(sys.argv[3], "w")
lflank = int(sys.argv[4])
rflank = int(sys.argv[5])


for line in inp_fp:
    k = line[0:line.find("\t")]
    footprint_vec = line[line.find("\t") + 1:] 
    
    per_orange = round(footprint_vec.count("F")*100/(len(footprint_vec) - 1), 3)
    real_lengths = footprint_length_dict[k]
    if len(real_lengths) == 0: # Naked DNA
        to_write = k + "\t" + str(per_orange) + "\t" + str(0)
        out_fp_per_orange_and_flen.write(to_write + "\n") 
    else:
  
        for l in real_lengths:
            to_write = k + "\t" + str(per_orange) + "\t" +  str(l)
            out_fp_per_orange_and_flen.write(to_write + "\n") 

out_fp_per_orange_and_flen.close() 
inp_fp.close() 
