import os
import sys
import re
from collections import defaultdict
from collections import OrderedDict
# it takes the input file from rule : plot_footprint_length_vs_per_orange 
# input: real_footprint_length_and_per_orange 
# sample file name: actual_footprint_length_in_clusters/merged_S2_to_reordered_20clust_147_50PNE_open_v1_peak_cluster_10_lf_14_rf_14_methylation_matrix_for_site_peak_5442_nclust_3_sam_flag_83~163_footprint_length.per_orange.tsv 

inp_fp = open(sys.argv[1])
footprint_fp = open(sys.argv[2]) 
out_fp = open(sys.argv[3], "w")

label_dict = {
    "Naked-DNA" : "0",
    "TF" : "1",
    "Nuc":  "2"
}

fpt_dict = OrderedDict()
for line in footprint_fp:
    k = line[0:line.find("\t")]
    footprint_str = line[line.find("\t") + 1: len(line) - 1]
    fpt_dict[k] = footprint_str 

labels_on_reads =  defaultdict(list)
for line in inp_fp:  
    line_length = len(line)
    d_loc = [m.start() for m in re.finditer("\t", line)]
    percent = float(line[d_loc[0] + 1: d_loc[1]] )
    flen = int(line[d_loc[1] + 1: line_length - 1])
    read_id = line[0:d_loc[0]]
    if percent <= 60: 
        labels_on_reads[read_id].append("Naked-DNA")
    elif percent>60 and flen < 100: 
        labels_on_reads[read_id].append("TF")
    else:
        labels_on_reads[read_id].append("Nuc")

for k in fpt_dict:
     if "TF" in labels_on_reads[k]: 
         to_write = k + "#" + label_dict["TF"] + "\t" + fpt_dict[k]
         out_fp.write(to_write + "\n")
     else:
         lab = labels_on_reads[k][0] # choose whatever comes first 
         to_write = k + "#" + label_dict[lab] + "\t" + fpt_dict[k]
         out_fp.write(to_write + "\n")
          
inp_fp.close()
out_fp.close() 
footprint_fp.close()
