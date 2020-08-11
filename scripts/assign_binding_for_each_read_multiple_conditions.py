import os
import sys
import re
import pickle
from collections import defaultdict
from collections import OrderedDict
import numpy as np 
from length_and_loc_with_absolute import get_real_footprint_length_with_abs_start
# it takes the input file from rule : plot_footprint_length_vs_per_orange 
# input: real_footprint_length_and_per_orange 
# sample file name: actual_footprint_length_in_clusters/merged_S2_to_reordered_20clust_147_50PNE_open_v1_peak_cluster_10_lf_14_rf_14_methylation_matrix_for_site_peak_5442_nclust_3_sam_flag_83~163_footprint_length.per_orange.tsv 

inp_fp = open(sys.argv[1])
footprint_fp = open(sys.argv[2]) 
out_fp_cond1 = open(sys.argv[3], "w")
out_fp_cond2 = open(sys.argv[4], "w")
out_fp_cond3 = open(sys.argv[5], "w")
out_fp_cond4 = open(sys.argv[6], "w")
out_fp_cond5 = open(sys.argv[7], "w")
out_fp_cond6 = open(sys.argv[8], "w")
per_methylation_vec_dict = pickle.load(open(sys.argv[9], "rb"))
out_fp_cond7 = open(sys.argv[10], "w")
out_fp_cond8 = open(sys.argv[11], "w")
occluded_dict = pickle.load(open(sys.argv[12], "rb"))
footprint_capped_dict = pickle.load(open(sys.argv[13], "rb"))

label_dict = {
    "Naked-DNA" : "0",
    "TF" : "1",
    "Nuc":  "2",
    "discard" : "3"
}

fpt_dict = OrderedDict()
fpt_abs_start_dict = defaultdict(list)
fpt_length_dict = defaultdict(list)
complete_fp_len_dict = defaultdict()
lflank =  int(footprint_fp.name.split("_lf_")[1].split("_rf_")[0])
rflank =  int(footprint_fp.name.split("_rf_")[1].split("_")[0]) 
for line in footprint_fp:
    d_loc = [m.start() for m in re.finditer("\t", line)]
    k = line[0:d_loc[0]]
    footprint_str = line[d_loc[0]+1:d_loc[1]]
    fpt_dict[k] = footprint_str 
    complete_footprint_vec = footprint_capped_dict[k]["footprint"] 
    complete_footprint_chr_start = int(footprint_capped_dict[k]["start"])
    complete_footprint_chr_end = int(footprint_capped_dict[k]["end"])
    fp_center = k.split(":")[1].split("-")[0]
    m_vec_start = int(fp_center)  - complete_footprint_chr_start -  lflank # absolute location on read 
    m_vec_stop = int(fp_center)  - complete_footprint_chr_start +  rflank + 1  # absolute location on read
    lvec, _, _, abs_loc = get_real_footprint_length_with_abs_start(footprint_str, m_vec_start, m_vec_stop, complete_footprint_vec)
    if "SRR3133329.32613656_32613656/1_overlapping" in k:
        print([footprint_str, m_vec_start, m_vec_stop, lvec, abs_loc, fp_center])
    fpt_abs_start_dict[k] = abs_loc 
    fpt_length_dict[k] = lvec
    complete_fp_len_dict[k] = len(complete_footprint_vec)
    
    
 

labels_on_reads_cond1 =  defaultdict(list)
labels_on_reads_cond2 =  defaultdict(list)
labels_on_reads_cond3 =  defaultdict(list)
labels_on_reads_cond4 =  defaultdict(list)
labels_on_reads_cond5 =  defaultdict(list)
labels_on_reads_cond6 =  defaultdict(list)
labels_on_reads_cond7 =  defaultdict(list)
labels_on_reads_cond8 =  defaultdict(list)
contribution_to_percentage_orange_cond8  =  defaultdict(list) 

for line in inp_fp:  
    line_length = len(line)
    d_loc = [m.start() for m in re.finditer("\t", line)]
    percent = float(line[d_loc[0] + 1: d_loc[1]] )
    flen = int(line[d_loc[1] + 1: d_loc[2]])
    read_id = line[0:d_loc[0]] 
    per_orange_for_fp = float(line[d_loc[2]+1: line_length - 1]) 
################################# condition 1 #####################
    if percent <= 60: 
        labels_on_reads_cond1[read_id].append("Naked-DNA")
    elif percent>60 and flen < 100: 
        labels_on_reads_cond1[read_id].append("TF")
    else:
        labels_on_reads_cond1[read_id].append("Nuc")

################################# condition 2 #####################
    if percent <= 60: 
        labels_on_reads_cond2[read_id].append("Naked-DNA")
    elif percent>60 and flen <=80: 
        labels_on_reads_cond2[read_id].append("TF")
    else:
        labels_on_reads_cond2[read_id].append("Nuc")
################################# condition 3 #####################
    if percent <= 60: 
        labels_on_reads_cond3[read_id].append("Naked-DNA")
    elif percent>60 and flen <=60: 
        labels_on_reads_cond3[read_id].append("TF")
    else:
        labels_on_reads_cond3[read_id].append("Nuc")
################################# condition 4 #####################
    if percent <= 60: 
        labels_on_reads_cond4[read_id].append("Naked-DNA")
    elif percent>60 and flen <=50: 
        labels_on_reads_cond4[read_id].append("TF")
    else:
        labels_on_reads_cond4[read_id].append("Nuc")
################################# condition 5 #####################
    if percent <= 60: 
        labels_on_reads_cond5[read_id].append("Naked-DNA")
    elif percent>60 and flen <=120: 
        labels_on_reads_cond5[read_id].append("TF")
    else:
        labels_on_reads_cond5[read_id].append("Nuc")

################################# condition 6 #####################
    if (percent <= 60) and float (per_methylation_vec_dict[read_id]["per_c"]) > 25: 
        labels_on_reads_cond6[read_id].append("Naked-DNA")
    elif (percent <= 60) and float (per_methylation_vec_dict[read_id]["per_c"]) <=25:
        labels_on_reads_cond6[read_id].append("Nuc")
   
    elif percent>60 and flen <=50: 
        labels_on_reads_cond6[read_id].append("TF")
    else:
        labels_on_reads_cond6[read_id].append("Nuc")

################################# condition 7 #####################
    if (percent <= 30) and float (per_methylation_vec_dict[read_id]["per_c"]) > 25: 
        labels_on_reads_cond7[read_id].append("Naked-DNA")
    elif (percent <= 30) and float (per_methylation_vec_dict[read_id]["per_c"]) <=25:
        labels_on_reads_cond7[read_id].append("Nuc")
   
    elif percent>30 and flen <=50: 
        labels_on_reads_cond7[read_id].append("TF")
    else:
        labels_on_reads_cond7[read_id].append("Nuc")

################################# condition 8 #####################
# use occluded to decide naked or occluded DNA 
# just add occluded condition and keep rest as condition 7
###################################################################
    max_of_edge = occluded_dict[read_id]["max_of_edge"]  
    pcap_total = occluded_dict[read_id]["percentage_capital"]
    #complete_footprint_vec = footprint_capped_dict[read_id]["footprint"] 
    #complete_footprint_chr_start = int(footprint_capped_dict[read_id]["start"])
    #complete_footprint_chr_end = int(footprint_capped_dict[read_id]["end"])
    
    # will have to check for edge condition
    # accept edged footprint only if they are greater than 50 to assign them a nucleomsome - they can't be assinged as TF because we don't have that information
    # 
    # get start and end of the 
    # chr3L:616610-616610^7`SRR3133326.2526820_2526820/1_overlapping`99~147	.....................	..Z..................	TACGATAATTGGTTTTTTTTT

    #dot_vec = "."*len(complete_footprint_vec)
    #if (complete_footprint_vec == dot_vec) and\
    #   (max_of_edge > 130) and (pcap_total !="NA"): 
    #   labels_on_reads_cond8[read_id].append("Nuc")
    if (percent <= 30) and float (per_methylation_vec_dict[read_id]["per_c"]) > 25: 
        labels_on_reads_cond8[read_id].append("Naked-DNA")
    elif (percent <= 30) and float (per_methylation_vec_dict[read_id]["per_c"]) <=25:
        labels_on_reads_cond8[read_id].append("Nuc")
   
    elif percent>30 and flen <=50 : 
        labels_on_reads_cond8[read_id].append("TF")
        contribution_to_percentage_orange_cond8[read_id].append(per_orange_for_fp)
    else:
        labels_on_reads_cond8[read_id].append("Nuc")
        contribution_to_percentage_orange_cond8[read_id].append(per_orange_for_fp)

for k in fpt_dict:
################################ condition 1 ##################################
    if "TF" in labels_on_reads_cond1[k]: 
        to_write = k + "#" + label_dict["TF"] + "\t" + fpt_dict[k]
        out_fp_cond1.write(to_write + "\n")
    else:
        lab = labels_on_reads_cond1[k][0] # choose whatever comes first 
        to_write = k + "#" + label_dict[lab] + "\t" + fpt_dict[k]
        out_fp_cond1.write(to_write + "\n")
############################### condition 2 ##################################
    if "TF" in labels_on_reads_cond2[k]:
        to_write = k + "#" + label_dict["TF"] + "\t" + fpt_dict[k]
        out_fp_cond2.write(to_write + "\n")      
    else: 
        lab = labels_on_reads_cond2[k][0] # choose whatever comes first 
        to_write = k + "#" + label_dict[lab] + "\t" + fpt_dict[k]
        out_fp_cond2.write(to_write + "\n")
############################### condition 3 ##################################
    if "TF" in labels_on_reads_cond3[k]:
        to_write = k + "#" + label_dict["TF"] + "\t" + fpt_dict[k]
        out_fp_cond3.write(to_write + "\n")      
    else: 
        lab = labels_on_reads_cond3[k][0] # choose whatever comes first 
        to_write = k + "#" + label_dict[lab] + "\t" + fpt_dict[k]
        out_fp_cond3.write(to_write + "\n")
############################### condition 4 ##################################
    if "TF" in labels_on_reads_cond4[k]:
        to_write = k + "#" + label_dict["TF"] + "\t" + fpt_dict[k]
        out_fp_cond4.write(to_write + "\n")      
    else: 
        lab = labels_on_reads_cond4[k][0] # choose whatever comes first 
        to_write = k + "#" + label_dict[lab] + "\t" + fpt_dict[k]
        out_fp_cond4.write(to_write + "\n")
############################### condition 5 ##################################
    if "TF" in labels_on_reads_cond5[k]:
        to_write = k + "#" + label_dict["TF"] + "\t" + fpt_dict[k]
        out_fp_cond5.write(to_write + "\n")      
    else: 
        lab = labels_on_reads_cond5[k][0] # choose whatever comes first 
        to_write = k + "#" + label_dict[lab] + "\t" + fpt_dict[k]
        out_fp_cond5.write(to_write + "\n")
############################### condition 6 ##################################
    if "TF" in labels_on_reads_cond6[k]:
        to_write = k + "#" + label_dict["TF"] + "\t" + fpt_dict[k]
        out_fp_cond6.write(to_write + "\n")      
    else: 
        lab = labels_on_reads_cond6[k][0] # choose whatever comes first 
        to_write = k + "#" + label_dict[lab] + "\t" + fpt_dict[k]
        out_fp_cond6.write(to_write + "\n")
############################### condition 7 ##################################
    if "TF" in labels_on_reads_cond7[k]:
        to_write = k + "#" + label_dict["TF"] + "\t" + fpt_dict[k]
        out_fp_cond7.write(to_write + "\n")      
    else: 
        lab = labels_on_reads_cond7[k][0] # choose whatever comes first 
        to_write = k + "#" + label_dict[lab] + "\t" + fpt_dict[k]
        out_fp_cond7.write(to_write + "\n")
############################### condition 8 ##################################
    #if "TF" in labels_on_reads_cond8[k]:
    #    to_write = k + "#" + label_dict["TF"] + "\t" + fpt_dict[k]
    #    out_fp_cond8.write(to_write + "\n")      
    #else: 
    #    lab = labels_on_reads_cond8[k][0] # choose whatever comes first 
    #    to_write = k + "#" + label_dict[lab] + "\t" + fpt_dict[k]
    #    out_fp_cond8.write(to_write + "\n")
    if len(contribution_to_percentage_orange_cond8[k]) == 0: # it should be naked
        to_write = k + "#" + label_dict["Naked-DNA"] + "\t" + fpt_dict[k]
        out_fp_cond8.write(to_write + "\n")
    else: 
        max_contributor_id = np.argmax(contribution_to_percentage_orange_cond8[k])
        abs_start_loc = fpt_abs_start_dict[k][max_contributor_id] 
        length_to_check = fpt_length_dict[k][max_contributor_id]
        if "SRR3133329.32613656_32613656/1_overlapping" in k:
            print ([abs_start_loc, length_to_check, complete_fp_len_dict[k]])
        lab = labels_on_reads_cond8[k][max_contributor_id] # choose whatever comes first 
        if (abs_start_loc == 0) and length_to_check <= 50:
            lab = "discard" # discard
        elif (abs_start_loc + length_to_check == complete_fp_len_dict[k]) and (length_to_check<=50):
            lab = "discard" # discard
        to_write = k + "#" + label_dict[lab] + "\t" + fpt_dict[k]
        out_fp_cond8.write(to_write + "\n")
        
         
          
inp_fp.close()
out_fp_cond1.close() 
out_fp_cond2.close() 
out_fp_cond3.close() 
out_fp_cond4.close() 
out_fp_cond5.close()
out_fp_cond6.close()
out_fp_cond7.close()
out_fp_cond8.close()
footprint_fp.close()
