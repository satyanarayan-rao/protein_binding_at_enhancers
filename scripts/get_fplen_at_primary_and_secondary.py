import os
import sys
from length_and_loc_with_absolute import get_real_footprint_length_with_abs_start
import pickle 
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

for line in inp_fp:
    # example line: 
    # chr3L:632728-632728^10`SRR3133326.1889480_1889480/1_overlapping`99~147	MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM...................................................................................................................................................................................................................................................................................................MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
    line_items = line.strip().split() 
    read_id = line_items[0]
    fp_str = line_items[1] 
    m_replaced_str = fp_str.replace("M", ".")
    primary_str_in_boundary = m_replaced_str[primary_lf: primary_rf]
    secondary_str_in_boundary = m_replaced_str[secondary_lf: secondary_rf]
    flen_list, _, loc_first, abs_start  = get_real_footprint_length_with_abs_start(
                   primary_str_in_boundary, primary_lf, primary_rf, 
                   m_replaced_str)
    #print(loc_first)
    max_of_edge = occluded_dict[read_id]["max_of_edge"]
    pcap_total = occluded_dict[read_id]["percentage_capital"]
    if (max_of_edge > 130) and (pcap_total != "NA"):
        to_write = "\t".join(["primary", str(max_of_edge), "both", read_id])
        out_fp.write(to_write + "\n") 
    elif len(flen_list) == 0:
        to_write = "\t".join(["primary", "0", "naked", read_id])
        out_fp.write(to_write + "\n") 
    else: 
        for l, abs_s in zip(flen_list, abs_start): 
            if abs_s + l > lextend + s_peak - lflank: 
                #print([abs_s, l, lextend, s_peak, lflank, read_id, "primary"])
                to_write = "\t".join(["primary", str(l), "both_fp", read_id]) 
                out_fp.write(to_write + "\n")
            else:
                to_write = "\t".join(["primary", str(l), "single", read_id]) 
                out_fp.write(to_write + "\n")
                

    flen_list, _, loc_first, abs_start  = get_real_footprint_length_with_abs_start(
                   secondary_str_in_boundary, secondary_lf, secondary_rf, 
                   m_replaced_str)
    if (max_of_edge > 130) and (pcap_total != "NA"):
        to_write = "\t".join(["secondary", str(max_of_edge), "both", read_id ])
        out_fp.write(to_write + "\n") 
    elif len(flen_list) == 0:
        to_write = "\t".join(["secondary", "0", "naked", read_id])
        out_fp.write(to_write + "\n") 
    else: 
        for l, abs_s in zip(flen_list, abs_start):
            #print (l + s_loc)
            if  abs_s  < lextend + rflank:
                #print([abs_s, l, lextend, s_peak, lflank, read_id, "secondary"])
                to_write = "\t".join(["secondary", str(l), "both_fp", read_id])
                out_fp.write(to_write + "\n")
            else:
                to_write = "\t".join(["secondary", str(l), "single", read_id])
                out_fp.write(to_write + "\n")
               
out_fp.close()
inp_fp.close()
