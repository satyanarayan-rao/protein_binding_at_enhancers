import os
import sys
import re
import pickle
from collections import defaultdict
from collections import OrderedDict
def get_real_footprint_length (m_vec, m_vec_start, m_vec_stop, complete_vec): 
    """
    extract complete footprint length from the selected methylation vector
    exmaple: 
    complete vec: . . . . F F F F . . . F F F . . . . . F F F F F 
    index       : 0 1 2 3 4 5 6 7 8 9 1011121314151617181920212223
    flank vec:  :             F F . . . F F F . . . . . F
                              6 7 8 9 10111213141516171819 
    output: [4, 3, 5]
    """
    # first get all footprint lengths in m_vec
    flen_lengths = [len(a) for a in m_vec.split(".")] 
    if m_vec.startswith("F"):
        flen_lengths[0] = flen_lengths[0] +  len(complete_vec[0:m_vec_start].split(".")[-1])
    if m_vec.endswith("F"):
        flen_lengths[-1] = flen_lengths[-1] +  len(complete_vec[m_vec_stop:].split(".")[0])
    # exclude all zeros
    return_list = [a for a in flen_lengths if a!=0] 
    # prepare a footprint length vector - at each index it will tell what is the footprint size that index it associated with
    first_f = False
    cnt = 0
    loc_first = []
    gap = True
    for c in m_vec:
        if c == 'F':
            first_f = True
            if gap == True:
                loc_first.append (cnt)
                gap = False
        else:
            first_f = False
            gap = True
        cnt +=1
    out_vec = [0]*len(m_vec)
    for i in range (len(loc_first)):
        for j in range(loc_first[i], len(m_vec)):
            if m_vec[j] == "F":
                out_vec[j] = return_list[i]
            else:
                break
    return return_list, out_vec, loc_first

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
    start = 0
    end = len(footprint_vec)
    len_vec, out_vec, start_loc = get_real_footprint_length(footprint_vec, start, end, complete_vec) 
    if len(len_vec) > 0:  
        for lv, sl in zip (len_vec, start_loc): 
            fstart = read_start + sl
            fend = read_start + lv + sl - 1 
            to_write = "\t".join ([chrom, str(fstart), str(fend), str(lv), read_id]) 
            out_fp.write(to_write + "\n") 
    else: # handling read without footprint
        to_write = "\t".join ([chrom, str(read_start), str(read_end), str(0), read_id]) 
        out_fp.write(to_write + "\n") 

out_fp.close()
inp_fp.close()
