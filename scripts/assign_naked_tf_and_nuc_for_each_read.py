import os
import sys
import re
from collections import defaultdict
# it takes the input file from rule : plot_footprint_length_vs_per_orange 
# input: real_footprint_length_and_per_orange 
# sample file name: actual_footprint_length_in_clusters/merged_S2_to_reordered_20clust_147_50PNE_open_v1_peak_cluster_10_lf_14_rf_14_methylation_matrix_for_site_peak_5442_nclust_3_sam_flag_83~163_footprint_length.per_orange.tsv 

inp_fp = open(sys.argv[1])
out_fp = open(sys.argv[2], "w")
labels_on_reads =  defaultdict(list)
out_read_and_label_fp = open(sys.argv[3], "w")
for line in inp_fp:  
    line_length = len(line)
    d_loc = [m.start() for m in re.finditer("\t", line)]
    percent = float(line[d_loc[0] + 1: d_loc[1]] )
    flen = int(line[d_loc[1] + 1: line_length - 1])
    read_id = line[0:d_loc[0]]
    if percent <= 60: 
        to_write = line[0: line_length - 1] + "\t" + "Naked-DNA"
        out_fp.write (to_write + "\n") 
        labels_on_reads[read_id].append("Naked-DNA")
    elif percent>60 and flen < 100: 
        to_write = line[0: line_length - 1] + "\t" + "TF" 
        labels_on_reads[read_id].append("TF")
        out_fp.write (to_write + "\n")
    else:
        to_write = line[0: line_length - 1] + "\t" + "Nuc" 
        out_fp.write (to_write + "\n")
        labels_on_reads[read_id].append("Nuc")

for k in labels_on_reads:
     if "TF" in labels_on_reads[k]: 
         to_write = k + "\t" + "TF" 
         out_read_and_label_fp.write(to_write + "\n")
     else:
         lab = labels_on_reads[k][0] # choose whatever comes first 
         to_write = k + "\t" + lab
         out_read_and_label_fp.write(to_write + "\n")
          
inp_fp.close()
out_fp.close() 
out_read_and_label_fp.close() 
