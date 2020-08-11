import os
import sys
import pickle
import re
from collections import defaultdict 
inp_fp = open(sys.argv[1]) 
occluded_dict = pickle.load(open(sys.argv[2], "rb")) 
out_fp = open(sys.argv[3], "w") 

read_flag_dict = defaultdict(lambda : False) 
for k in occluded_dict:
    read_flag_dict[k] = True

cnt = 0 
for line in inp_fp:
    # example line: head -1 extended_binding_labeled_reads_for_peak/suppressed_merged_S2_to_reordered_20clust_147_50PNE_open_v1_peak_cluster_9_lf_15_rf_15_site_peak_229_sam_flag_99~147_extend_from_peak_left_150_right_150.cond8.tsv
    # chr2L:480305-480305^9`SRR3133326.2671172_2671172/1_overlapping`99~147#1	..............................................................................................................................FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF........MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
    # example from dict 
    # python $NGS_SCRIPTS_DIR/view_sample_from_pkl.py occluded_edges_on_methylation_vec/suppressed_merged_S2_to_all_open_and_closed_mnase_peaks_lf_15_rf_15_min_flen_10_wobble_gap_1_occluded.pkl 
    # ['chr3R:19975740-19975740^8`SRR3133326.2817484_2817484/1_overlapping`83~163', {'total_vec_len': 278, 'total_letters_on_read': 30, 'percentage_capital': 0.53, 'max_of_edge': 133, 'start': 145, 'end': 277}]

    if cnt % 3 == 0:  # only footprint entry changes
        tab_loc = line.find("\t")
        read_id = line[0:line.find("\t")].split("#")[0]
        fp_vec = line[tab_loc+1:len(line) - 1]
        if read_flag_dict[read_id] == True:
            if (occluded_dict[read_id]['max_of_edge'] > 130) and\
               (occluded_dict[read_id]['percentage_capital'] !="NA") and\
               (fp_vec.count('F') == 0):  # used e only when legnth of edge is greater than 130 
                # find the last M from left 
                idx = 0 
                for c in fp_vec: 
                    if c!='M':
                        break
                    else:
                        idx +=1 
                e_start = idx + occluded_dict[read_id]["start"]
                e_end = idx + occluded_dict[read_id]["end"] + 1
                e_vec = fp_vec[0:idx] + fp_vec[idx:e_start] + "E"*(e_end - e_start) + fp_vec[e_end: len(fp_vec)] 
                e_vec_301 = e_vec[0:301] 
                #out_fp.write(str(idx) + "\t" + str (occluded_dict[read_id]["start"]) + "\t" + str(occluded_dict[read_id]["end"]) + "\t" + str(e_start) + "\t" + str(e_end) + "\t" +  line[0:tab_loc] + "\t" + e_vec + "\n") 
                out_fp.write(line[0:tab_loc] + "\t" + e_vec_301 + "\n") 
            else:
                out_fp.write(line)
        else:
            out_fp.write(line)
    else:
        out_fp.write(line)
    cnt +=1 


out_fp.close()
inp_fp.close()
