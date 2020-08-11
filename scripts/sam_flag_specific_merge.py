import sys 
from collections import defaultdict

sam_99_fp = open(sys.argv[1])
sam_83_fp = open(sys.argv[2])
out_fp = open(sys.argv[3], "w")
sam_flag_wise_dict = defaultdict(list)
tmp_83_dict = defaultdict(list)
all_binding_labels = []
for line in sam_99_fp:
    # example line: head -1 extended_binding_labeled_reads_for_peak/suppressed_merged_S2_to_reordered_20clust_147_50PNE_open_v1_peak_cluster_9_lf_15_rf_15_site_peak_229_sam_flag_99~147_extend_from_peak_left_150_right_150.cond8.footprint_ordered.tsv
    # chr2L:480305-480305^9`SRR3133327.34990438_34990438/1_overlapping`99~147#0	MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM......................................................................FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
    binding_label = int(line.split("\t")[0].split("#")[-1])
    sam_flag_wise_dict[binding_label].append(line) 
    if len(all_binding_labels) == 0: 
        all_binding_labels.append(binding_label) 
    elif binding_label not in all_binding_labels:
        all_binding_labels.append(binding_label)
    else: 
        continue

for line in sam_83_fp:
    # example line: head -1 extended_binding_labeled_reads_for_peak/suppressed_merged_S2_to_reordered_20clust_147_50PNE_open_v1_peak_cluster_9_lf_15_rf_15_site_peak_229_sam_flag_83~163_extend_from_peak_left_150_right_150.cond8.footprint_ordered.tsv
    # chr2L:480305-480305^9`SRR3133326.2818_2818/1_overlapping`83~163#0	MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM.......................................................................................................................................................................
    binding_label = int(line.split("\t")[0].split("#")[-1])
    tmp_83_dict[binding_label].append(line) 
    if len(all_binding_labels) == 0: 
        all_binding_labels.append(binding_label) 
    elif binding_label not in all_binding_labels:
        all_binding_labels.append(binding_label)
    else: 
        continue
# sort the 83 flag
for k in tmp_83_dict:
    #rev_list = tmp_83_dict[k][::-1]
    rev_list = tmp_83_dict[k]
    for j in rev_list:
        sam_flag_wise_dict[k].append(j)


sorted_binding_labels = sorted(all_binding_labels)
for k in sorted_binding_labels:
    for l in sam_flag_wise_dict[k]:
        out_fp.write(l) 

out_fp.close()
sam_99_fp.close()
sam_83_fp.close()       

