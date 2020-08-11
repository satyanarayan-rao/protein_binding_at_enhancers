import os
import sys
import re
from collections import defaultdict 
label_dict = {
    "0" : "Naked-DNA",
    "1" : "TF",
    "2" : "Nuc", 
    "3" : "discard"
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
            total_reads +=1 
            binding_id = f_line[f_line.find("#")+1: f_line.find("\t")] 
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


header = "file_name_without_sam_flag\ttotal_reads\ttotal_naked\ttotal_tf\ttotal_nuc\toccup_naked\toccup_tf\toccup_nuc\ttotal_discard\toccup_discard"
out_fp.write (header + "\n") 

count_dict_83 = get_counts_for_each_flag (inp_fp_83, sam_flag = "83~163")
count_dict_99 = get_counts_for_each_flag (inp_fp_99, sam_flag = "99~147")
for k in count_dict_83: 
     fname=k
     total = count_dict_83[k]["total"]  + count_dict_99[k]["total"]
     total_naked = count_dict_83[k]["Naked-DNA"]  + count_dict_99[k]["Naked-DNA"]
     total_tf = count_dict_83[k]["TF"]  + count_dict_99[k]["TF"]
     total_nuc = count_dict_83[k]["Nuc"]  + count_dict_99[k]["Nuc"]
     per_naked = round ((total_naked/total)*100, 2) 
     per_tf = round ((total_tf/total)*100, 2)
     per_nuc = round ((total_nuc/total)*100, 2) 
     total_discard = count_dict_83[k]["discard"] + count_dict_99[k]["discard"]  
     per_discard = round((total_discard/total)*100, 2) 
     to_write = "\t".join([fname, 
                str(total), str(total_naked), str(total_tf), str(total_nuc),
                str(per_naked), str(per_tf), str(per_nuc), str(total_discard), str(per_discard)])
     out_fp.write (to_write + "\n") 
     
out_fp.close()
