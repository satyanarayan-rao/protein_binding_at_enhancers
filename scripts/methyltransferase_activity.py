import os
import sys
from collections import defaultdict
import re
inp_fp = open(sys.argv[1])
context_dict = defaultdict(lambda : 0)
total_cnt = 0
for line in inp_fp:
    bs_seq = line[0:line.find("\t")]
    mvec = line[line.find("\t") +1: len(line) - 1 ] 
    # find all capital letters - as in the ones got methylated 
    m_loc = []
    cnt = -1
    for c in mvec: 
        cnt +=1 
        if c == ".":
            continue
        elif c == c.upper():
            m_loc.append(cnt) 
        else:
            continue
    for v in m_loc:
        if v == 0:
            if mvec[0] == "Z":
                context_dict["CG"] +=1
                total_cnt +=1
            else:
                context_dict["at_start"] +=1
                total_cnt +=1
        else:
            if mvec[v] == "Z": 
                context_dict["CG"] +=1 
                total_cnt +=1
            elif bs_seq[v-1:v+1] == "GC":
                context_dict["GC"] +=1 
                total_cnt +=1
            else:
                context_dict["HCH"] +=1 
                total_cnt +=1
              

inp_fp.close()

for k in ["CG", "GC", "HCH", "at_start"]:
    print ("%5d\t%0.3f\t%5s"%(context_dict[k], (context_dict[k]/total_cnt)*100, k))


