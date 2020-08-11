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
    accepted_labels = list(map(str, range(9)))
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
            if binding_id in accepted_labels:
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


header = "file_name_without_sam_flag\ttotal_tf_left\ttotal_tf_right\ttotal_tf_tf\ttotalreads\tobs_by_exp\tdistance_between_peaks" 

out_fp.write (header + "\n") 

count_dict_83 = get_counts_for_each_flag (inp_fp_83, sam_flag = "83~163")
count_dict_99 = get_counts_for_each_flag (inp_fp_99, sam_flag = "99~147")

for k in count_dict_83: 
     total_list = [] 
     percent_list = [] 
     fname=k
     total = count_dict_83[k]["total"]  + count_dict_99[k]["total"]
     tf_tf = count_dict_83[k]["occup_T_T"] + count_dict_99[k]["occup_T_T"]
     total_in_left = 0
     total_in_right = 0 
     dist_val = k.split("_sec_peak_at_")[1].split("_")[0]
     for j in [3,4,5]: # left tf 
         total_in_left += count_dict_83[k][label_dict[str(j)]] +  count_dict_99[k][label_dict[str(j)]] 
     for j in [1,4,7]: # right_tf
         total_in_right += count_dict_83[k][label_dict[str(j)]] +  count_dict_99[k][label_dict[str(j)]]

     obs_by_exp = "NA" 
     if (total_in_left >0) and (total_in_right >0):
        expected = (total_in_left/total)*(total_in_right/total)  
        obs_by_exp = round((tf_tf/total)/expected,3) 
     to_write = "\t".join([k, str(total_in_left), 
                           str(total_in_right), str(tf_tf),
                           str(total), str(obs_by_exp), dist_val])
     out_fp.write (to_write + "\n") 
     
out_fp.close()
