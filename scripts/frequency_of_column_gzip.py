import os
import sys
import re
import gzip 
from collections import defaultdict

inp_fp = gzip.open (sys.argv[1])
out_fp = open (sys.argv[2], "w")
column_id = int(sys.argv[3])
max_lim = int(sys.argv[4])
ftype = sys.argv[5]
bed = sys.argv[6]
cnt_dict = defaultdict (lambda : 0 ) 

total_cnt = 0
tab_gz = re.compile(bytes("\t", encoding = "ascii"))
for line in inp_fp:
     d_loc = [m.start () for m in re.finditer(tab_gz, line)]
     flen = bytes.decode(line[d_loc[column_id - 2] + 1: d_loc[column_id - 1]])
     cnt_dict[flen] += 1 
     total_cnt += 1 

for i in range (max_lim + 1):
    frac = round(cnt_dict[str(i)]/total_cnt, 5)
    to_write = "\t".join([str(i), str(cnt_dict[str(i)]), str (frac), ftype, bed]) 
    out_fp.write(to_write + "\n")
frac_na = round(cnt_dict["NA"]/total_cnt, 5)
na_write = "\t".join(["NA", str(cnt_dict["NA"]), str (frac_na), ftype, bed]) 
out_fp.write(na_write + "\n")
out_fp.close()
inp_fp.close()
