import os
import sys
from collections import defaultdict
import re

inp_fp = open(sys.argv[1])
cnt_dict = defaultdict(lambda : defaultdict(lambda : 0))
is_read_counted_for_peak = defaultdict(lambda : defaultdict(lambda : False))
for line in inp_fp:
    # example line
    # chr3R	19748037	19748134	peak_3545	chr3R:19748119-19748119^3	-	0	67	1	chr3R	19747966	19748249	SRR3133326.1459798_1459798/1_overlapping`83~163	.
    d_loc = [m.start() for m in re.finditer("\t", line)]
    k = "\t".join([line[d_loc[7]+1:d_loc[8]],
                   line[d_loc[2]+1:d_loc[3]], 
                   line[d_loc[5]+1:d_loc[7]],
                   line[d_loc[3]+1:d_loc[5]]]) 
    read_id = line[d_loc[-2] +1: d_loc[-1]] 
    if ("83~163" in line) and (is_read_counted_for_peak[k][read_id] == False): 
        cnt_dict[k]["83~163"] +=1
        is_read_counted_for_peak[k][read_id] = True
    elif ("99~147" in line) and (is_read_counted_for_peak[k][read_id] == False):
        cnt_dict[k]["99~147"] +=1
        is_read_counted_for_peak[k][read_id] = True
        
    
for key in cnt_dict: 
    print ("{}\t{}\t{}".format(key, cnt_dict[key]["99~147"], cnt_dict[key]["83~163"])) 

inp_fp.close()
