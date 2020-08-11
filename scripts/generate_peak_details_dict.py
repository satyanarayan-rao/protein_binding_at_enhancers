import os
import sys
from collections import defaultdict
import pickle
input_fp = open (sys.argv[1]) 
out_fp = open(sys.argv[2], "wb")

out_dict = defaultdict(dict)
header = input_fp.readline().strip().split()
total = len(header)
for line in input_fp:
    line_items = line.strip().split() 
    for idx in range(1, total):
        out_dict[line_items[0]][header[idx]] = line_items[idx] 

pickle.dump(dict(out_dict), out_fp)
out_fp.close()
