import os
import sys
from collections import defaultdict
import pickle
import re

inp_fp = open(sys.argv[1])
len_dict = pickle.load(open(sys.argv[2], "rb"))
out_fp = open(sys.argv[3], "w")

for line in inp_fp:
    k = line[0:line.find("#")]
    f_str = line[line.find("\t") + 1: len(line) - 1] 
    per_o = round ((f_str.count ("F"))*100/len(f_str), 3)
    if per_o == 0:
        to_write = line[0: line.find("\t")] + "\t" + "0" + "\t" + "0"
        out_fp.write(to_write + "\n")
    else:
        for v in len_dict[k]:
            to_write = line[0:line.find("\t")] + "\t" + str(per_o) + "\t" + str(v)
            out_fp.write(to_write + "\n")

out_fp.close()
inp_fp.close()
