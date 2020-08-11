import os 
import sys
from collections import defaultdict
import re
bed_fp = open(sys.argv[1]) 
out_fp = open(sys.argv[2], "w")

# chr2L	11806729	11806730	chr2L:11806729-11806729^20	.	+	chr2L	11806644	11806923	SRR3133326.203452_203452/1_overlapping`83~163	.	.h..h.H..x...h.H.Z.H.........Z....x...x........X...........h......h..zxhhhh.........z.....................xh...h.....hh..x...x..xhh.......h..h.h........h....h..z..........h.................X......Z.H.......h.....Zx....H.............h........Z....Z..........z.....X................
peak_read_cnt_dict = defaultdict(lambda : defaultdict(lambda : 0 ))
for line in bed_fp:
    d_loc = [m.start() for m in re.finditer("\t", line)]
    chr_loc_and_cl = line[d_loc[2] + 1: d_loc[3]] 
    sam_flag = line[d_loc[8] + 1: d_loc[9]].split("`")[-1]  
    peak_read_cnt_dict[chr_loc_and_cl][sam_flag]  += 1 


for k in peak_read_cnt_dict:
    mnase_cl = k.split ("^")[-1]
    chr_loc = k.split ("^")[0]
    to_write = "%40s\t%5s\t%5s\t%5s"%(chr_loc, peak_read_cnt_dict[k]["99~147"],
                             peak_read_cnt_dict[k]["83~163"], mnase_cl)
    out_fp.write(to_write + "\n")
