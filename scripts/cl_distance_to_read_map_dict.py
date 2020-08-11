import os
import sys
import re
import pickle
from collections import defaultdict
inp_fp = open(sys.argv[1])
out_fp = open(sys.argv[2], "wb")

distance_to_read_map_dict = defaultdict(lambda: defaultdict(list))
for line in inp_fp:
    # example line
    # chr3R	19748037	19748134	peak_3545	chr3R:19748119-19748119^3	-	0	67	1	chr3R	19747966	19748249	SRR3133326.1459798_1459798/1_overlapping`83~163	.
    d_loc = [m.start() for m in re.finditer("\t", line)] 
    # prepare key using `peak_id`, `secondary_peak_loc`, `strand` and `sam_flag` 
    k = "\t".join([
             line[d_loc[2]+1:d_loc[3]],
             line[d_loc[6]+1:d_loc[7]],
             line[d_loc[4]+1:d_loc[5]]]) 
    if "83~163" in line:
        distance_to_read_map_dict[k]["83~163"].append(line[d_loc[11]+1:d_loc[12]])
    elif "99~147" in line:
        distance_to_read_map_dict[k]["99~147"].append(line[d_loc[11]+1:d_loc[12]]) 
    
pickle.dump(dict(distance_to_read_map_dict), out_fp)
out_fp.close()
inp_fp.close()
    
    

