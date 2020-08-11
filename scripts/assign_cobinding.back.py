import os
import sys
from length_and_loc_with_absolute import get_real_footprint_length_with_abs_start
import pickle 
from collections import defaultdict
inp_fp = open(sys.argv[1])
lflank = int(sys.argv[2])
rflank = int(sys.argv[3])
s_peak = int(sys.argv[4])
lextend = int(sys.argv[5])

primary_lf = 0 + lextend - lflank  
primary_rf = 0 + lextend + rflank + 1

secondary_lf = 0 + lextend + s_peak - lflank
secondary_rf = 0 + lextend + s_peak + rflank + 1
out_fp = open(sys.argv[6], "w")
occluded_dict = pickle.load(open(sys.argv[7],"rb"))
out_fp_150 = open(sys.argv[8], "w")
primay_assignment_dict = defaultdict(dict)
secondary_assignment_dict = defaultdict(dict)
ordered_reads_on_cobinding = defaultdict(lambda : defaultdict(list))

pair_map_dict = {
    "0-0": 0, 
    "0-1": 1, 
    "0-2": 2, 
    "1-0": 3, 
    "1-1": 4, 
    "1-2": 5, 
    "2-0": 6, 
    "2-1": 7, 
    "2-2": 8 
}

for line in inp_fp:
    # example line: 
    # chr3L:632728-632728^10`SRR3133326.1889480_1889480/1_overlapping`99~147	MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM...................................................................................................................................................................................................................................................................................................MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
    line_items = line.strip().split() 
    read_id = line_items[0]
    fp_str = line_items[1] 
    m_replaced_str = fp_str.replace("M", ".")
    primary_str_in_boundary = m_replaced_str[primary_lf: primary_rf]
    secondary_str_in_boundary = m_replaced_str[secondary_lf: secondary_rf]
    primary_per_orange = round (primary_str_in_boundary.count ('F')/len (primary_str_in_boundary), 3)*100 
    secondary_per_orange = round (secondary_str_in_boundary.count ('F')/len (secondary_str_in_boundary), 3)*100 
    flen_list, _, loc_first, abs_start  = get_real_footprint_length_with_abs_start(
                   primary_str_in_boundary, primary_lf, primary_rf, 
                   m_replaced_str)
    max_of_edge = occluded_dict[read_id]["max_of_edge"]
    pcap_total = occluded_dict[read_id]["percentage_capital"]
    if (max_of_edge > 130) and (pcap_total != "NA") and primary_per_orange<=30:
        primay_assignment_dict[read_id]["binding_label"] = "2" # nuc
    elif len(flen_list) == 0:
        primay_assignment_dict[read_id]["binding_label"] = "0" # naked
    elif primary_per_orange <=30:
        primay_assignment_dict[read_id]["binding_label"] = "0" # naked
       
    elif primary_per_orange >30: 
        is_lt_fplen = False
        for l in flen_list:
            if l <=100:
                is_lt_fplen = True
        if is_lt_fplen == True: 
            primay_assignment_dict[read_id]["binding_label"] = "1"
        else:
            primay_assignment_dict[read_id]["binding_label"] = "2"
           
    flen_list, _, loc_first, abs_start  = get_real_footprint_length_with_abs_start(
                   secondary_str_in_boundary, secondary_lf, secondary_rf, 
                   m_replaced_str)
    if (max_of_edge > 130) and (pcap_total != "NA") and secondary_per_orange<=30:
        secondary_assignment_dict[read_id]["binding_label"] = "2" # nuc
    elif len(flen_list) == 0:
        secondary_assignment_dict[read_id]["binding_label"] = "0" # naked
    elif secondary_per_orange <=30:
        secondary_assignment_dict[read_id]["binding_label"] = "0" # naked
    elif secondary_per_orange>30: 
        is_lt_fplen = False
        for l in flen_list:
            if l <=100:
                is_lt_fplen = True
        if is_lt_fplen == True: 
            secondary_assignment_dict[read_id]["binding_label"] = "1"
        else:
            secondary_assignment_dict[read_id]["binding_label"] = "2"
    one_of_nine_labels = pair_map_dict[primay_assignment_dict[read_id]["binding_label"] + "-" + secondary_assignment_dict[read_id]["binding_label"]] 

    ordered_reads_on_cobinding[one_of_nine_labels]["read_id"].append(read_id)
    ordered_reads_on_cobinding[one_of_nine_labels]["vec"].append(fp_str)

all_labels = list(ordered_reads_on_cobinding.keys())
all_labels.sort()


for k in all_labels:
    for read_id, fp_str in zip(ordered_reads_on_cobinding[k]["read_id"], 
        ordered_reads_on_cobinding[k]["vec"]):
        to_write = read_id + "#"  + str(k) + "\t" + fp_str
        fp_150 = fp_str[int((len(fp_str) - 1)/2) - 150:int((len(fp_str) - 1)/2) + 150 + 1]
        to_write_150 = read_id + "#"  + str(k) + "\t" + fp_150
        out_fp.write(to_write + "\n")
        out_fp_150.write(to_write_150 + "\n") 
    
out_fp.close()
inp_fp.close()
