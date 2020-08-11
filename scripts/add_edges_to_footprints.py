import os
import sys
import pickle
import re
from collections import defaultdict 
inp_fp = open(sys.argv[1]) 
out_fp = open(sys.argv[2], "w") 

def get_border_methylated_cytosines (m_vec, len_vec):
    left_border = None
    right_border = None
    cnt = 0
    for c in m_vec: 
        if (c==".") or (c == c.lower()):
            cnt+=1 
        elif c == c.upper():
            break 
    left_border = cnt # point at the capital letter in python index scheme
    cnt = 0
    for c in m_vec[::-1]:
        if (c==".") or (c == c.lower()):
            cnt+=1 
        elif c == c.upper():
            break         
    right_border = len(m_vec) - cnt - 1 # point at the capital letter in python index scheme
    
    if left_border == len_vec:
        left_border = None
        right_border = None
    return left_border, right_border    
     

#read_flag_dict = defaultdict(lambda : False) 
#for k in occluded_dict:
#    read_flag_dict[k] = True

cnt = 0 
for line in inp_fp:
    # example line: head -1 suppressed_merged_S2_to_reordered_20clust_147_50PNE_open_v1_peak_lf_15_rf_15_tf_footprint_min_length_10_wobble_gap_1.bed
    #chr2L	11806729	11806730	chr2L:11806729-11806729^20	.	+	chr2L	11806460	11806747	SRR3133326.952549_952549/1_overlapping`99~147	.	..........................................................................................FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF........................................................	..X..Z.....Z.............H..............................................H...Z.......H....H............................z....z...x....x....z....z.....z....h.....................................x........z...H.......z...................H...................Z......H.......z....................	GGCTGCGTTTGCGTGGAAAGAAGAGCAAAAATATTTTAATTTTAATTGGAAAGGGAATGGGAATTGGAATTGCTTTCGGGGGGGCATTGCAAATTTTATATTTTTAATTATTTTAATTTGAAATGAGTAGTGTTGTGTGGTGTGGTTTTGGTGTTAATAAATCAAATTTATTGTATTTATATTTAGAAGAGTTTTATGTGTGTGCATTTTTTTGATTAGTTTGTTTTTTTTGCTTTTTTTCTAGATTTTAGACGGGGGGCCTTTTTTTGTTTTTTTTTTTTTATTATT 
    # example from dict 
    # python view_pkl.py suppressed_merged_S2_to_reordered_20clust_147_50PNE_open_v1_peak_lf_15_rf_15_min_flen_10_wobble_gap_1_occluded.pkl 
    # ['chr2L:11806729-11806729^20`SRR3133326.952549_952549/1_overlapping`99~147', {'total_vec_len': 288, 'total_letters_on_read': 24, 'percentage_capital': 0.5, 'max_of_edge': 28, 'start': 260, 'end': 287}]
    d_loc = [m.start() for m in re.finditer("\t", line)]  
    

    fp_vec = line[d_loc[10]+1:d_loc[11]]
    m_vec = line[d_loc[11]+1:d_loc[12]]  
    l_mvec = len(m_vec) 
    left_border, right_border = get_border_methylated_cytosines(m_vec=m_vec, len_vec=l_mvec)
    new_vec = None
    if (left_border !=None) and (right_border!=None):
        new_vec = "F"*(left_border) + fp_vec[left_border: right_border+1] + 'F'*(l_mvec - right_border - 1) 
        #print ([len (fp_vec) - len(new_vec), fp_vec, new_vec, m_vec]) 
        #print(len(fp_vec) - len(new_vec))
    else:
        #No capital letters found - just assign the whole read as footprint 
        new_vec =  'F'*(l_mvec) 
        
    to_write = "\t".join([line[0:d_loc[10]], new_vec, line[d_loc[11]+1:]])
    out_fp.write(to_write)

out_fp.close()
inp_fp.close()
