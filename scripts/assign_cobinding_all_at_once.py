import os
import sys
import pickle 
import numpy as np
from collections import defaultdict
from modules_for_cobinding import count_m_at_start_end
from modules_for_cobinding import add_edge_colors
from modules_for_cobinding import get_per_orange_for_each_footprint   
from modules_for_cobinding import get_read_start_and_end_from_lex_rex_reads 
from modules_for_cobinding import get_count_and_percentage_methylation
from length_and_loc_with_absolute import get_real_footprint_length_with_abs_start
from length_and_loc import get_real_footprint_length  

#-------------------------------------# 
occluded_dict = pickle.load(open(sys.argv[1], "rb"))
job_id_pkl = pickle.load(open(sys.argv[2], "rb")) 
job_counter = 0
for job_k in job_id_pkl:

    inp_fp = open(job_id_pkl[job_k][0])
    lflank = int(job_id_pkl[job_k][1])
    rflank = int(job_id_pkl[job_k][2])
    s_peak = int(job_id_pkl[job_k][3])
    lextend = int(job_id_pkl[job_k][4])
    out_fp = open(job_id_pkl[job_k][5], "w")
    #occluded_dict = pickle.load(open(job_id_pkl[job_k][7],"rb"))
    out_fp_150 = open(job_id_pkl[job_k][7], "w")
    cobinding_verbose_fp = open(job_id_pkl[job_k][8])
    mvec_out_fp = open(job_id_pkl[job_k][9], "w")
    mvec_out_150 = open(job_id_pkl[job_k][10], "w")
    primary_assignment_dict = defaultdict(dict)
    secondary_assignment_dict = defaultdict(dict)
 
    #inp_fp = open(sys.argv[1])
    #lflank = int(sys.argv[2])
    #rflank = int(sys.argv[3])
    #s_peak = int(sys.argv[4])
    #lextend = int(sys.argv[5])
    #out_fp = open(sys.argv[6], "w")
    #occluded_dict = pickle.load(open(sys.argv[7],"rb"))
    #out_fp_150 = open(sys.argv[8], "w")
    
    primary_lf = 0 + lextend - lflank  
    primary_rf = 0 + lextend + rflank + 1
    
    secondary_lf = 0 + lextend + s_peak - lflank
    secondary_rf = 0 + lextend + s_peak + rflank + 1
    primay_assignment_dict = defaultdict(dict)
    secondary_assignment_dict = defaultdict(dict)
    ordered_reads_on_cobinding = defaultdict(lambda : defaultdict(list))
    to_color_edges = defaultdict(dict)
    
    cobinding_dict = defaultdict(lambda : False)
    pair_map_dict = {
        "0-0": 0, 
        "0-1": 1, 
        "0-2": 2, 
        "1-0": 3, 
        "1-1": 4, 
        "1-2": 5, 
        "2-0": 6, 
        "2-1": 7, 
        "2-2": 8,
        "3-0": 9, # discard on left naked on right
        "3-1": 10, # discard on left TF on right
        "3-2": 11,
        "0-3": 12, 
        "1-3": 13,
        "2-3": 14,
        "3-3": 15   
    }
    ############## create cobinding footprint mvec and bsseq dict ############ 
    cobinding_verb_dict = defaultdict(dict)
    cnt = 0 
    line_mod_to_type_map_dict = {
        0: "footprint",
        1: "mvec",
        2: "bs_seq" }
    for line in cobinding_verbose_fp: 
        tab_loc = line.find("\t")
        cobinding_verb_dict[line[0:tab_loc]][line_mod_to_type_map_dict[cnt%3]] = line[tab_loc:len(line) - 1] 
        cnt +=1
    
    
    ##################

    
    for line in inp_fp:
        # example line: 
        # head -1 filtered_cobinding_reads/suppressed_merged_S2_to_all_open_and_closed_mnase_peaks_cluster_199_lf_15_rf_15_sec_peak_at_55_lex_300_rex_300_site_peak_229_3_sam_flag_83~163.tsv
        # chr2L:480305-480305^199`SRR3133326.2247651_2247651/1_overlapping`83~163	MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMFFFFFFF..........................................................................FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF...................................................FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF.FFFFFFFFFMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
        line_items = line.strip().split() 
        read_id = line_items[0]
        fp_str = line_items[1] 
        m_counts_start, m_counts_end = count_m_at_start_end(fp_str) 
        read_start_abs, read_end_abs = get_read_start_and_end_from_lex_rex_reads(fp_str) 
        m_replaced_str = fp_str.replace("M", ".")
        primary_str_in_boundary = m_replaced_str[primary_lf: primary_rf]
        secondary_str_in_boundary = m_replaced_str[secondary_lf: secondary_rf]
        primary_per_orange = round (primary_str_in_boundary.count ('F')/len (primary_str_in_boundary), 3)*100 
        secondary_per_orange = round (secondary_str_in_boundary.count ('F')/len (secondary_str_in_boundary), 3)*100 
        flen_list_primary, _, loc_first_primary, abs_start_primary  = get_real_footprint_length_with_abs_start(
                       primary_str_in_boundary, primary_lf, primary_rf, 
                       m_replaced_str) 
        flen_list_secondary, _, loc_first_secondary, abs_start_secondary  = get_real_footprint_length_with_abs_start(
                       secondary_str_in_boundary, secondary_lf, secondary_rf, 
                       m_replaced_str) 
        one_of_nine_labels = None
        
        ########### Get percentage orange per footprint for each peak ########### 
        # will be used in making decision on what read to consider ##############  
        
        primary_percentage_orange_per_footprint = get_per_orange_for_each_footprint(primary_str_in_boundary)
        secondary_percentage_orange_per_footprint = get_per_orange_for_each_footprint(secondary_str_in_boundary) 
          
        # get percentage methylation in the range +15 to -30 and -15 to +30 and choose the minimum ############ 
        primary_left_neighbor_mvec = cobinding_verb_dict[read_id]["mvec"][primary_lf - 15:primary_rf] 
        primary_right_neighbor_mvec = cobinding_verb_dict[read_id]["mvec"][primary_lf:primary_rf + 15]  
        secondary_left_neighbor_mvec = cobinding_verb_dict[read_id]["mvec"][secondary_lf - 15:secondary_rf] 
        secondary_right_neighbor_mvec = cobinding_verb_dict[read_id]["mvec"][secondary_lf:secondary_rf + 15]  
          
        primary_per_c_l, primary_total_upper_l, primary_total_lower_l, primary_total_l = get_count_and_percentage_methylation (primary_left_neighbor_mvec) 
        primary_per_c_r, primary_total_upper_r, primary_total_lower_r, primary_total_r = get_count_and_percentage_methylation (primary_right_neighbor_mvec) 
        
        secondary_per_c_l, secondary_total_upper_l, secondary_total_lower_l, secondary_total_l = get_count_and_percentage_methylation (secondary_left_neighbor_mvec) 
        secondary_per_c_r, secondary_total_upper_r, secondary_total_lower_r, secondary_total_r = get_count_and_percentage_methylation (secondary_right_neighbor_mvec) 
        #print ("@@@@@@@@@@@@@@@@@@@@@@")
        #print ([primary_str_in_boundary, secondary_str_in_boundary, primary_percentage_orange_per_footprint, secondary_percentage_orange_per_footprint ])
        ######################################
    
        ################# start assigning the label using the data collected above ###
        
        ## first check for cobinding to be assined to long read as it is special in a way that primary and secondary alone conditions don't apply to this ###### 
        # start with % orange winned concept at both peaks   ############  
        if (len(primary_percentage_orange_per_footprint) > 0 ) and\
           (len(secondary_percentage_orange_per_footprint) > 0): 
            # covers TF-TF : single long footprint
            # covers TF-TF : independent footprints making TF-TF  
            # covers NUC-NUC: independent/same footprints making NUC-NUC
            # covers TF-Naked: based on nonzero percentage 
            # covers Naked-TF: based on nonzero percentage 
            # covers NUC-Naked: based on nonzero percentage
            # covers Naked-NUC: based on nonzero percentage
            # covers Naked-Naked: based on nonzero percentage
            max_orange_of_primary = np.argmax(primary_percentage_orange_per_footprint)
            max_orange_primary = primary_percentage_orange_per_footprint [max_orange_of_primary] 
            length_of_max_orange_read_primary = flen_list_primary[max_orange_of_primary]  
            abs_start_of_max_orange_read_primary = abs_start_primary[max_orange_of_primary]  
    
            max_orange_of_secondary = np.argmax(secondary_percentage_orange_per_footprint)  
            max_orange_secondary = secondary_percentage_orange_per_footprint [max_orange_of_secondary] 
            length_of_max_orange_read_secondary = flen_list_secondary[max_orange_of_secondary]  
            abs_start_of_max_orange_read_secondary = abs_start_secondary[max_orange_of_secondary]          
            # case when it is actually a long footprint representing  # 
            # both percentge orange primary and secondary >30% ########
            # cobound footprint length <=100 ########################## 
            # cobound footprint is not edge  ########################## 
            # abs start obtained from both primary and sec to be same #
            #print (["here", abs_start_of_max_orange_read_primary, abs_start_of_max_orange_read_secondary, length_of_max_orange_read_primary, length_of_max_orange_read_secondary, read_start_abs, read_end_abs, max_orange_primary, max_orange_secondary])
            if (max_orange_primary > 30) and (max_orange_secondary > 30):
                if (abs_start_of_max_orange_read_primary == abs_start_of_max_orange_read_secondary):
                    if (abs_start_of_max_orange_read_primary != read_start_abs) and\
                       (abs_start_of_max_orange_read_primary + length_of_max_orange_read_primary - 1!= read_end_abs):
                        if s_peak <=100: # d = 70, 15 from site1 and 15 from site2 
                            if length_of_max_orange_read_primary <=100:
                                cobinding_dict[read_id] = True
                                primary_assignment_dict[read_id]["binding_label"] = "1" # TF
                                secondary_assignment_dict[read_id]["binding_label"] = "1" # TF 
                            else:
                                primary_assignment_dict[read_id]["binding_label"] = "2" # Nuc
                                secondary_assignment_dict[read_id]["binding_label"] = "2" # Nuc 
                        else: # here it will come only if footprint length is greater than 100 because abs start is same and sec peak is (d+30) >100 apart : d>70
                            primary_assignment_dict[read_id]["binding_label"]  = "2" # NUC 
                            secondary_assignment_dict[read_id]["binding_label"]  = "2" # NUC
    
                    elif length_of_max_orange_read_primary <=50:
                         primary_assignment_dict[read_id]["binding_label"] = "3" # Discard 
                         secondary_assignment_dict[read_id]["binding_label"] = "3" # Discard
                    else:
                         primary_assignment_dict[read_id]["binding_label"] = "2" # Nuc
                         secondary_assignment_dict[read_id]["binding_label"] = "2" # Nuc
                             
                else: # two footprints are different handle, them separately - primary first
                    if (abs_start_of_max_orange_read_primary != read_start_abs) and\
                       (abs_start_of_max_orange_read_primary + length_of_max_orange_read_primary - 1!= read_end_abs): 
                        if length_of_max_orange_read_primary <=50:
                            primary_assignment_dict[read_id]["binding_label"] =  "1"  # TF - as it is not the edge
                        else:
                            primary_assignment_dict[read_id]["binding_label"] =  "2"  # NUC 
                    else: # its an edge 
                        if length_of_max_orange_read_primary >50: 
                            primary_assignment_dict[read_id]["binding_label"] = "2" # NUC
                        else: 
                            primary_assignment_dict[read_id]["binding_label"] = "3" # Discard
                               
                    
                    if (abs_start_of_max_orange_read_secondary != read_start_abs) and\
                       (abs_start_of_max_orange_read_secondary + length_of_max_orange_read_secondary - 1!= read_end_abs): 
                        if length_of_max_orange_read_secondary <=50:
                            secondary_assignment_dict[read_id]["binding_label"] =  "1"  # TF
                        else:
                            secondary_assignment_dict[read_id]["binding_label"] =  "2"  # NUC
                    else:
                        if length_of_max_orange_read_secondary > 50: 
                            secondary_assignment_dict[read_id]["binding_label"] = "2" # NUC
                        else:
                            secondary_assignment_dict[read_id]["binding_label"] = "3" # Discard 
                            
            elif (max_orange_primary >30) and (max_orange_secondary <= 30):
                # check for the edges always 
                # first work on secondary
                if (abs_start_of_max_orange_read_secondary != read_start_abs) or\
                   (abs_start_of_max_orange_read_secondary + length_of_max_orange_read_secondary -1 != read_end_abs): 
                    secondary_assignment_dict[read_id]["binding_label"] = "0" # naked 
                else:
                    secondary_assignment_dict[read_id]["binding_label"] = "3" # discard
                # note: here I can't assign NUC for secondary because percentage is <=30
                # handle the primary condition now : in this case it can only be TF, NUC or discard
                if (abs_start_of_max_orange_read_primary == read_start_abs) or\
                   (abs_start_of_max_orange_read_primary + length_of_max_orange_read_primary - 1 == read_end_abs):
                    if length_of_max_orange_read_primary > 50: 
                        primary_assignment_dict[read_id]["binding_label"] = "2" # NUC
                    else:
                        primary_assignment_dict[read_id]["binding_label"] = "3" # discard
                else: # genunine footprint : within read
                    if length_of_max_orange_read_primary <=50: 
                        primary_assignment_dict[read_id]["binding_label"] = "1" # TF
                    else:
                        primary_assignment_dict[read_id]["binding_label"] = "2" # NUC
            elif (max_orange_primary<=30) and (max_orange_secondary > 30):
                if (abs_start_of_max_orange_read_primary != read_start_abs) or\
                   (abs_start_of_max_orange_read_primary + length_of_max_orange_read_primary -1 != read_end_abs): 
                    primary_assignment_dict[read_id]["binding_label"] = "0" # naked 
                else:
                    primary_assignment_dict[read_id]["binding_label"] = "3" # discard
                # note: here I can't assign NUC for secondary because percentage is <=30
                # handle the primary condition now : in this case it can only be TF, NUC or discard
                if (abs_start_of_max_orange_read_secondary == read_start_abs) or\
                   (abs_start_of_max_orange_read_secondary + length_of_max_orange_read_secondary - 1 == read_end_abs):
                    if length_of_max_orange_read_secondary > 50: 
                        secondary_assignment_dict[read_id]["binding_label"] = "2" # NUC
                    else:
                        secondary_assignment_dict[read_id]["binding_label"] = "3" # discard
                else: # genunine footprint : within read
                    if length_of_max_orange_read_secondary <=50: 
                        secondary_assignment_dict[read_id]["binding_label"] = "1" # TF
                    else:
                        secondary_assignment_dict[read_id]["binding_label"] = "2" # NUC               
            else: 
                if (abs_start_of_max_orange_read_primary != read_start_abs) or\
                   (abs_start_of_max_orange_read_primary + length_of_max_orange_read_primary -1 != read_end_abs): 
                    primary_assignment_dict[read_id]["binding_label"] = "0" # naked 
                else:
                    primary_assignment_dict[read_id]["binding_label"] = "3" # discard
                # note: here I can't assign NUC for secondary because percentage is <=30
                # handle the primary condition now : in this case it can only be TF, NUC or discard
                if (abs_start_of_max_orange_read_secondary != read_start_abs) or\
                   (abs_start_of_max_orange_read_secondary + length_of_max_orange_read_secondary - 1 != read_end_abs):
                                       
                    secondary_assignment_dict[read_id]["binding_label"] = "0"
                else:
                    secondary_assignment_dict[read_id]["binding_label"] = "3"
                
        elif (len(primary_percentage_orange_per_footprint) > 0 ) and\
             (len(secondary_percentage_orange_per_footprint) == 0):
            # footprint found in primary but no footprint in secondary
            # covers TF-Naked : 
            # covers Naked-Naked:
            # covers Nuc-Naked: 
            # covers discard-Naked: 
            
            # coming to think of it - after fixing the edges I think no footprint in preak +- will always orginate from naked case - I guess no need to check for percetntage methylation in the neighbor
            # directly assign naked to secondary with the above logic
            secondary_assignment_dict[read_id]["binding_label"] = "0" 
            # handle the logic for primary
            max_orange_of_primary = np.argmax(primary_percentage_orange_per_footprint)
            max_orange_primary = primary_percentage_orange_per_footprint [max_orange_of_primary] 
            length_of_max_orange_read_primary = flen_list_primary[max_orange_of_primary]  
            abs_start_of_max_orange_read_primary = abs_start_primary[max_orange_of_primary]  
            # check for the edge as always
            if max_orange_primary>30:
                 # check for TF, NUC, edge
                if (abs_start_of_max_orange_read_primary == read_start_abs) or\
                   (abs_start_of_max_orange_read_primary + length_of_max_orange_read_primary - 1 == read_end_abs):
                    if length_of_max_orange_read_primary > 50: 
                        primary_assignment_dict[read_id]["binding_label"] = "2" # NUC
                    else:
                        primary_assignment_dict[read_id]["binding_label"] = "3" # discard
                elif length_of_max_orange_read_primary > 50:
                        primary_assignment_dict[read_id]["binding_label"] = "2" # NUC
                else:
                        primary_assignment_dict[read_id]["binding_label"] = "1" # TF
                  
     
            else:
                 primary_assignment_dict[read_id]["binding_label"] = "0"
        elif (len(primary_percentage_orange_per_footprint) == 0 ) and\
             (len(secondary_percentage_orange_per_footprint) > 0):
            # footprint found in primary but no footprint in secondary
            # covers Naked-TF : 
            # covers Naked-Naked:
            # covers Naked-Nuc: 
            # covers Naked-discard: 
            
            # coming to think of it - after fixing the edges I think no footprint in preak +- will always orginate from naked case - I guess no need to check for percetntage methylation in the neighbor
            primary_assignment_dict[read_id]["binding_label"] = "0" 
            # handle the logic for secondary now
            max_orange_of_secondary = np.argmax(secondary_percentage_orange_per_footprint)
            max_orange_secondary = secondary_percentage_orange_per_footprint [max_orange_of_secondary] 
            length_of_max_orange_read_secondary = flen_list_secondary[max_orange_of_secondary]  
            abs_start_of_max_orange_read_secondary = abs_start_secondary[max_orange_of_secondary]  
            # check for the edge as always
            if max_orange_secondary>30:
                 # check for TF, NUC, edge
                if (abs_start_of_max_orange_read_secondary == read_start_abs) or\
                   (abs_start_of_max_orange_read_secondary + length_of_max_orange_read_secondary - 1 == read_end_abs):
                    if length_of_max_orange_read_secondary > 50: 
                        secondary_assignment_dict[read_id]["binding_label"] = "2" # NUC
                    else:
                        secondary_assignment_dict[read_id]["binding_label"] = "3" # discard
                elif length_of_max_orange_read_secondary > 50:
                        secondary_assignment_dict[read_id]["binding_label"] = "2" # NUC
                else:
                        secondary_assignment_dict[read_id]["binding_label"] = "1" # TF
                  
     
            else:
                 secondary_assignment_dict[read_id]["binding_label"] = "0"        
        else:    
            primary_assignment_dict[read_id]["binding_label"] = "0"        
            secondary_assignment_dict[read_id]["binding_label"] = "0"        
        #print ([read_id, primary_assignment_dict[read_id],secondary_assignment_dict[read_id] ]) 
        #print([inp_fp.name, read_id])
        label_comb = primary_assignment_dict[read_id]["binding_label"] + "-" +\
                     secondary_assignment_dict[read_id]["binding_label"] 
        one_of_nine_labels =  pair_map_dict[label_comb]
    
        ordered_reads_on_cobinding[one_of_nine_labels]["read_id"].append(read_id)
        ordered_reads_on_cobinding[one_of_nine_labels]["vec"].append(fp_str)
        ordered_reads_on_cobinding[one_of_nine_labels]["mvec"].append(cobinding_verb_dict[read_id]["mvec"])
        
    
    all_labels = list(ordered_reads_on_cobinding.keys())
    all_labels.sort()
    
    
    for k in all_labels:
        for read_id, fp_str, mvec_str in zip(
            ordered_reads_on_cobinding[k]["read_id"], 
            ordered_reads_on_cobinding[k]["vec"],
            ordered_reads_on_cobinding[k]["mvec"]):
            modified_fp_str = fp_str
            to_write = read_id + "#"  + str(k) + "\t" + modified_fp_str
            fp_150 = modified_fp_str[int((len(modified_fp_str) - 1)/2) - 150:int((len(modified_fp_str) - 1)/2) + 150 + 1]
            to_write_150 = read_id + "#"  + str(k) + "\t" + fp_150
            out_fp.write(to_write + "\n")
            out_fp_150.write(to_write_150 + "\n") 
            
    
            to_write = read_id + "#"  + str(k) + "\t" + mvec_str
            mvec_out_fp.write(to_write + "\n") 
            mvec_str_150 = mvec_str[int((len(modified_fp_str) - 1)/2) - 150:int((len(modified_fp_str) - 1)/2) + 150 + 1]
            to_write = read_id + "#"  + str(k) + "\t" + mvec_str_150
            mvec_out_150.write(to_write + "\n")
            
        
    out_fp.close()
    inp_fp.close()
    mvec_out_fp.close()
    mvec_out_150.close()
    print ("Finished jon {} ...".format(job_k))
    #if job_counter == 3:
    #   break
    #job_counter +=1 
