import os
import sys
import re
from collections import defaultdict 
label_dict = {
    "0" : "occup_D_D",
    "1" : "occup_D_T",
    "2" : "occup_D_N", 
    "3" : "occup_T_D",
    "4" : "occup_T_T",
    "5" : "occup_T_N",
    "6" : "occup_N_D",
    "7" : "occup_N_T",
    "8" : "occup_N_N"
}
accepted_binding_states = list(map(str, range(9)))
def get_counts_for_each_flag (flist_name, sam_flag=None):
    """
    flist_name: filename with 83 or 99 flag
    return: a dictionary with:
          key: actual tsv file name without sam flag
              sub_key: <0, 1, 2> for naked, TF, Nuc
                value: <count_0, count_1, count_2>  
                value: total
    """
    fp_sam = open(flist_name)
    count_dict = {} 
    to_replace_flag = "sam_flag_" + sam_flag  + "_"
    for line in fp_sam:
        fname = line.strip() 
        # example fname
        # binding_labeled_reads_for_peak/suppressed_merged_S2_to_reordered_20clust_147_50PNE_open_v1_peak_cluster_6_lf_15_rf_15_site_peak_2057_sam_flag_83~163_binding_label.cond4.tsv
        fp = open(fname)
        total_reads = 0 
        occup_count_dict = defaultdict(lambda : 0) 
        for f_line in fp:
            binding_id = f_line[f_line.find("#")+1: f_line.find("\t")] 
            if binding_id in accepted_binding_states: 
                total_reads +=1 
                occup_count_dict[label_dict[binding_id]] +=1 
        occup_count_dict["total"] = total_reads
        fp.close()
        fname_as_key = fname.replace(to_replace_flag, "")
        count_dict[fname_as_key] = occup_count_dict     
    fp_sam.close()    
    return count_dict
     
inp_fp_83 = sys.argv[1]
inp_fp_99 = sys.argv[2]
out_fp = open(sys.argv[3], "w")


header_fname_and_total = "file_name_without_sam_flag\ttotal_reads"

header_states_count = "\t".join(["total_" + label_dict[str(i)] for i in range(len(label_dict)) ])
header_states_percent = "\t".join(["per_" + label_dict[str(i)] for i in range(len(label_dict)) ])
header = "\t".join ([header_fname_and_total, header_states_count, header_states_percent]) 

out_fp.write (header + "\n") 

count_dict_83 = get_counts_for_each_flag (inp_fp_83, sam_flag = "83~163")
count_dict_99 = get_counts_for_each_flag (inp_fp_99, sam_flag = "99~147")

for k in count_dict_83: 
     total_list = [] 
     percent_list = [] 
     fname=k
     total = count_dict_83[k]["total"]  + count_dict_99[k]["total"]
     for j in range(len(label_dict)): 
         total_in_state = count_dict_83[k][label_dict[str(j)]] +  count_dict_99[k][label_dict[str(j)]] 
         percent_in_state = round((total_in_state/total)*100,2) 
         total_list.append(total_in_state)
         percent_list.append(percent_in_state) 
     total_state_cnt_str = "\t".join(map(str,total_list)) 
     total_state_per_str = "\t".join(map(str,percent_list)) 
     
     to_write = "\t".join([fname, 
                str(total), total_state_cnt_str, total_state_per_str])
     out_fp.write (to_write + "\n") 
     
out_fp.close()
