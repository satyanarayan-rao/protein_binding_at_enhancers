import os
import sys
import re

inp_fp = open(sys.argv[1])
cutoff = int(sys.argv[2])
out_fp = open(sys.argv[3], "w")

for line in inp_fp:
    d_loc = [m.start() for m in re.finditer("\t", line)]
    footprint_str = line[d_loc[10] + 1:d_loc[11]] 
    after_cutoff = []
    cnt_flag = 0
    cnt = 0
    for c in footprint_str:
        if c == "F":
            cnt_flag +=1
            after_cutoff.append(c)
        else:
            if cnt_flag >=cutoff:
                after_cutoff.append(c)
                cnt_flag = 0
            elif cnt_flag == 0:
                after_cutoff.append(c)
            else:
                for i in range(cnt - cnt_flag - 1, cnt):
                    after_cutoff[i] = "."
                cnt_flag = 0
                after_cutoff.append(c)
        cnt += 1
    if (cnt_flag !=0) and cnt_flag <cutoff:
        for i in range(len(footprint_str) - cnt_flag - 1, len(footprint_str)):
            after_cutoff[i] = "." 
    to_write = line [0:d_loc[10]] + "\t" + "".join(after_cutoff) + "\t" + line[d_loc[11] + 1:]
    out_fp.write(to_write)
    
inp_fp.close() 
out_fp.close()
