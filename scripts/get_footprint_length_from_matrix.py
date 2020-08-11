import os
import sys
import re

input_fp = open(sys.argv[1]) 
out_fp = open(sys.argv[2], "w") 

for line in input_fp:
    footprint_str = line[line.find("\t") + 1: len(line) - 1]
    cnt = 0
    for c in footprint_str:
        if c == "F": 
            cnt += 1
        else:
            if cnt > 0: 
                to_write = str(cnt)
                out_fp.write(to_write + "\n")
                cnt  = 0
    if cnt > 0: 
        to_write = str(cnt)
        out_fp.write(to_write + "\n")
out_fp.close()
 
    
