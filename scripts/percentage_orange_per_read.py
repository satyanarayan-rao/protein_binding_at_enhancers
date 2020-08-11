import os
import sys
import re

input_fp = open(sys.argv[1]) 
out_fp = open(sys.argv[2], "w") 

for line in input_fp:
    footprint_str = line[line.find("\t") + 1: len(line) - 1] 
    footprint_id = line[0:line.find("\t")]
    per_orange = round(footprint_str.count("F")*100/(len(footprint_str)), 3)
    out_fp.write(footprint_id + "\t" + str(per_orange) + "\n")

out_fp.close()
