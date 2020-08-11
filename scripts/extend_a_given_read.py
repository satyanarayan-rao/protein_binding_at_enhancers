import os
import sys
import pickle
import re
# head -1 footprint.bed
#chr2L	11806729	11806730	chr2L:11806729-11806729^20	.	+	chr2L	11806460	11806747	SRR3133326.952549_952549/1_overlapping`99~147`99~147	.	........F...................FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF......................FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF.....FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF.....FFFFFFFFFFFFFFFFFFFFFFF.....FFF.....FFFFFFF......................................

# head -1 merged_S2_to_reordered_20clust_147_50PNE_open_v1_peak_cluster_7_lf_22_rf_19_methylation_matrix_for_site_peak_3484_nclust_3_sam_flag_83~163.tsv 
# chr3R:5416975-5416975^7`SRR3133326.727267_727267/1_overlapping`83~163`83~163#4	...................FFFFFFFFFFFFFF.....FFFF

def extend_footprint_from_center (line, footprint_dict,
        lflank, rflank, lextend, rextend, strand): 
    key_without_hash = line[0:line.find("\t")]
    existing_footprint = line[line.find("\t") + 1: len(line) - 1] 
    complete_footprint = footprint_dict[key_without_hash]["footprint"]
    complete_mvec = footprint_dict[key_without_hash]["mvec"]
    complete_bsseq = footprint_dict[key_without_hash]["bs_seq"]
    complete_footprint_start = int(footprint_dict[key_without_hash]["start"])
    complete_footprint_end = int(footprint_dict[key_without_hash]["end"])
    peak_center = int(key_without_hash.split("-")[0].split(":")[-1])
    m_vec_start = peak_center - complete_footprint_start -  lflank
    m_vec_stop = peak_center - complete_footprint_start  + rflank + 1
    methylation_string = complete_footprint[m_vec_start: m_vec_stop]
    ####
    # lf = -3, rf = 3 
    # m_vec_start = 26
    #  m_vec_stop = 33 
    # extend_start = m_vec_start - (lextend - lf) 
    # extend_stop = m_vec_stop + (rextend - rf)  
    # 
    # complete_footprint_start and complete_footprint_end are in closed form length of footprint = 288, `11806747 - 11806460 + 1` 
    # actual  2324252627282930313233343536 
    #          . . . . . F F F F . . . F F 
    #               -3-2-1 0 1 2 3
    # desired -5 7 
    #          . . . . . F F F F . . . F F
    #           -5-4-3-2-1 0 1 2 3 4 5 6 7 
    # 
    #### 
    #lflank = int(sys.argv[3])
    #rflank = int(sys.argv[4])
    #lextend = int(sys.argv[5])
    #rextend = int(sys.argv[6])
    #    if strand == "-": 
    #        lexetend = int(sys.argv[6])
    #        rextend = int(sys.argv[5])
    #        lflank = int(sys.argv[4])
    #        rflank = int(sys.argv[3]) 
        
    desired_length = rextend - (0 - lextend) + 1 
    extend_start = peak_center - (lextend) # close interval
    extend_stop = peak_center + (rextend + 1 )  # open 
    extended_footprint = ""
    extended_mvec = ""
    extended_bsseq = ""
    #print ([extend_start, extend_stop, complete_footprint_start, complete_footprint_end])
    #sys.exit(-1)
    if (extend_start >= complete_footprint_start) and  (extend_stop <= complete_footprint_end + 1): 
        #print ([key_with_hash, "condition 1"])
        if strand == "+":
            extended_footprint = complete_footprint[extend_start - complete_footprint_start: extend_stop - complete_footprint_start]
            extended_mvec = complete_mvec[extend_start - complete_footprint_start: extend_stop - complete_footprint_start]
            extended_bsseq = complete_bsseq[extend_start - complete_footprint_start: extend_stop - complete_footprint_start]
        elif strand == "-": 
            extended_footprint = complete_footprint[extend_start - complete_footprint_start: extend_stop - complete_footprint_start][::-1]
            extended_mvec = complete_mvec[extend_start - complete_footprint_start: extend_stop - complete_footprint_start][::-1]
            extended_bsseq = complete_bsseq[extend_start - complete_footprint_start: extend_stop - complete_footprint_start][::-1]
            # recall that complete_footprint is on the watson strand always
            #print(["in", "in", key_without_hash,  extended_footprint, len(extended_footprint)])
    elif (extend_start >= complete_footprint_start) and  (extend_stop > complete_footprint_end + 1): 
        #print ([key_with_hash, "condition 2"])
        if strand == "+":
            extended_footprint = complete_footprint[extend_start - complete_footprint_start: complete_footprint_end - complete_footprint_start + 1]
            extended_mvec = complete_mvec[extend_start - complete_footprint_start: complete_footprint_end - complete_footprint_start + 1]
            extended_bsseq = complete_bsseq[extend_start - complete_footprint_start: complete_footprint_end - complete_footprint_start + 1]
            # add `M` of how much short in the end 
            extended_footprint = extended_footprint + "M"*(desired_length - len(extended_footprint))
            extended_mvec = extended_mvec + "M"*(desired_length - len(extended_mvec))
            extended_bsseq = extended_bsseq + "M"*(desired_length - len(extended_bsseq))
        elif strand == "-":
            extended_footprint = complete_footprint[extend_start - complete_footprint_start: complete_footprint_end - complete_footprint_start + 1][::-1] 
            extended_mvec = complete_mvec[extend_start - complete_footprint_start: complete_footprint_end - complete_footprint_start + 1][::-1] 
            extended_bsseq = complete_bsseq[extend_start - complete_footprint_start: complete_footprint_end - complete_footprint_start + 1][::-1] 
            #print([peak_center, extend_start, complete_footprint_start, complete_footprint_end])
            extended_footprint = "M"*(desired_length - len(extended_footprint)) + extended_footprint
            extended_mvec = "M"*(desired_length - len(extended_mvec)) + extended_mvec 
            extended_bsseq = "M"*(desired_length - len(extended_bsseq)) + extended_bsseq 
            #print(["in", "out", key_without_hash,  extended_footprint, len(extended_footprint)])
    elif (extend_start < complete_footprint_start) and (extend_stop <= complete_footprint_end + 1): 
        if strand == "+":
            extended_footprint = complete_footprint[0: extend_stop - complete_footprint_start]
            extended_mvec = complete_mvec[0: extend_stop - complete_footprint_start]
            extended_bsseq = complete_bsseq[0: extend_stop - complete_footprint_start]
            extended_footprint = "M"*(desired_length - len(extended_footprint)) + extended_footprint 
            extended_mvec = "M"*(desired_length - len(extended_mvec)) + extended_mvec 
            extended_bsseq = "M"*(desired_length - len(extended_bsseq)) + extended_bsseq
        #print(["out", "in", key_without_hash,  extended_footprint, len(extended_footprint)]) 
        elif strand == "-":
            extended_footprint = complete_footprint[0: extend_stop - complete_footprint_start][::-1]
            extended_mvec = complete_mvec[0: extend_stop - complete_footprint_start][::-1]
            extended_bsseq = complete_bsseq[0: extend_stop - complete_footprint_start][::-1]
            extended_footprint = extended_footprint + "M"*(desired_length - len(extended_footprint)) 
            extended_mvec = extended_mvec + "M"*(desired_length - len(extended_mvec)) 
            extended_bsseq = extended_bsseq + "M"*(desired_length - len(extended_bsseq)) 
            #print(["out", "in", key_without_hash,  extended_footprint, len(extended_footprint)]) 
        #print ([key_with_hash, "condition 3", extended_mvec])
        
    elif (extend_start < complete_footprint_start) and (extend_stop > complete_footprint_end + 1): 
        #print ([key_with_hash, "condition 4"])
        if strand == "+":
            to_add_in_left = "M"*(complete_footprint_start - extend_start)
            to_add_in_right = "M"*(extend_stop - (complete_footprint_end + 1) )
            extended_footprint = to_add_in_left + complete_footprint + to_add_in_right   
            extended_mvec = to_add_in_left + complete_mvec + to_add_in_right   
            extended_bsseq = to_add_in_left + complete_bsseq + to_add_in_right   
        elif strand == "-":
            to_add_in_left = "M"*(complete_footprint_start - extend_start)
            to_add_in_right = "M"*(extend_stop - (complete_footprint_end + 1) )
            extended_footprint = to_add_in_right + complete_footprint[::-1] + to_add_in_left
            extended_mvec = to_add_in_right + complete_mvec[::-1] + to_add_in_left
            extended_bsseq = to_add_in_right + complete_bsseq[::-1] + to_add_in_left
            #print(["out", "out", key_without_hash, extended_footprint, len(extended_footprint)])
    
    #print ([extend_start, extend_stop, complete_footprint_start, complete_footprint_end, extended_footprint])
    return extended_footprint, extended_mvec, extended_bsseq
