# python  scripts/report2pkl.py <*report.txt.gz> <*report.pkl>
import os
import sys
import pickle
import gzip
import re

input_fp = gzip.open (sys.argv[1])
dict_fp = open (sys.argv[2], "wb")
methylation_level_dict = {} 
for line in input_fp:
    line_items = bytes.decode(line)
    line_items = line_items.strip().split()
    dict_key = "`".join(line_items[0:2])
    methylated = int (line_items[3]) 
    unmethylated = int (line_items[4]) 
    percnt = round (float(methylated/(unmethylated+methylated)) * 100, 4)
    methylation_level_dict[dict_key] = percnt

pickle.dump(methylation_level_dict, dict_fp)
input_fp.close() 
dict_fp.close()
