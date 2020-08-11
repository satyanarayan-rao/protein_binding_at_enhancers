import sys
import re 
from collections import defaultdict 

inp_fp = open(sys.argv[1])
th = int(sys.argv[2])
#out_fp = open(sys.argv[3], "w")
peak_id_counter_dict =  defaultdict(lambda : 0)
peak_id_loc_dict =  defaultdict(list)
for line in inp_fp:
    # example line: head -1 all_open_enh_read_cnt_per_peak.tsv_ge_15.tsv
    # chr3R:19975666-19975666	   47	  115	peak_4746
    l_items = line.strip().split() 
    peak_id = l_items[-1] 
    peak_id_counter_dict[peak_id] +=1 
    peak_iter_id = peak_id + "_" + str(peak_id_counter_dict[peak_id])
    loc = int(l_items[0].split("-")[-1])
    if len(peak_id_loc_dict[peak_id]) == 0:
        peak_id_loc_dict[peak_id].append(loc)
        to_write = "\t".join([l_items[0], l_items[1], l_items[2], peak_iter_id ])
        print (to_write) 
    else:
        is_lt_th = False
        for l in peak_id_loc_dict[peak_id]: 
            if abs(loc - l) <= th:
                is_lt_th = True
        if is_lt_th == False:
            peak_id_loc_dict[peak_id].append(loc)
            to_write = "\t".join([l_items[0], l_items[1], l_items[2], peak_iter_id ])
            print (to_write) 
        else:
            continue
            
