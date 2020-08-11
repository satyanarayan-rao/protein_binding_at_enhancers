import os
import sys
import re
import pickle
from collections import defaultdict
from collections import OrderedDict
from length_and_loc import get_real_footprint_length
from length_and_loc_with_absolute import get_real_footprint_length_with_abs_start 
from capital_peracentage_and_unbound_stretches  import capital_percentage_and_stretch_of_unbound

inp_fp = open(sys.argv[1]) 
out_fp = open(sys.argv[2], "w")

for line in inp_fp:
    # example line:
    # chr3L	632728	632729	chr3L:632728-632728^10	.	+	chr3L	632480	632750	SRR3133326.1892464_1892464/1_overlapping`83~163	.	.................................................................FFFFFFFFFFF...................................................................................................................................................................................................	....................................H.Z.H.Z.....................Z.....z.....X....h...........................................................zx.z...z.z......z........z.........z..............h....h..........zx.........h......................h.............................	AATATATATACATATAATACATTCATAAATTTCCTTGCGTGCGTCATAATTTTCTATAATTCTCGATCTCATTACTGCCTAACCATTTCCATATCTACTTACAATATATACATATATATATATATATATTCCTAATCTATCAACAATCATATTAATCATACTCATCAAATTTCCCCAAAAATAACTAAATAACATAACTACTTTCTCAACTATTCTAAACTTATTTATATATAAACTTAAAACTAATCATAAATTAATTCTAATTCATAAA
    d_loc = [m.start() for m in re.finditer("\t", line)]  
    read_start = int(line[d_loc[6] +1: d_loc[7]]) 
    read_end = int(line[d_loc[7] +1: d_loc[8]]) 
    chrom = line[d_loc[5] +1: d_loc[6]]
    read_id = line[d_loc[8] +1: d_loc[9]]
    footprint_vec = line[d_loc[10] + 1:d_loc[11]]
    complete_vec = footprint_vec
    mvec = line[d_loc[11] + 1: d_loc[12]]
    start = 0
    end = len(footprint_vec)
    len_vec, _, loc_first, abs_loc = get_real_footprint_length_with_abs_start(footprint_vec, 0, end, complete_vec) 
    
    if len(len_vec) > 0:  
        for lv, sl in zip (len_vec, abs_loc): 
            if (sl == 0) or (lv + sl) == end: 
                fstart = read_start + sl - 1
                fend = read_start + lv + sl - 1 
                to_write = "\t".join ([chrom, str(fstart), str(fend), str(lv), read_id, "edge"]) 
                out_fp.write(to_write + "\n") 
            else:
                fstart = read_start + sl - 1
                fend = read_start + lv + sl - 1 
                to_write = "\t".join ([chrom, str(fstart), str(fend), str(lv), read_id, "inside"]) 
                out_fp.write(to_write + "\n")              
    else: # handling read without footprint
        to_write = "\t".join ([chrom, str(read_start), str(read_end), str(0), read_id, "zero"]) 
        out_fp.write(to_write + "\n") 

out_fp.close()
inp_fp.close()
