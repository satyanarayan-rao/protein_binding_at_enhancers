import os
import sys
import re
from collections import defaultdict
import pickle
from extend_a_given_read import extend_footprint_from_center
from capital_peracentage_and_unbound_stretches import capital_percentage_and_stretch_of_unbound
# Look at the rule get_real_footprint_lengths_and_oraange_on_read
def get_count_and_percentage_methylation (m_vec):
    total = 0
    total_lower = 0
    total_upper = 0
    per_c = 0
    for c in m_vec:
        if c == ".":
            continue
        elif c == c.lower(): 
            total +=1
            total_lower +=1
        elif c == c.upper():
            total +=1 
            total_upper +=1 
    if total>0:
        per_c = str(round(total_upper/total, 3)*100)
    else:
        per_c = "NA"
    return per_c, total_upper, total_lower, total

############# get %orange for each footprint ##############  
def get_per_orange_for_each_footprint (fvec):
    flen_list = [] 
    cnt = 0 
    f_switch = False
    fvec_len = len(fvec)
    for c in fvec:
        if c == 'F':
            cnt +=1
            f_switch = True
        else:
            if cnt!=0:
            	flen_list.append(round((cnt/fvec_len)*100, 2))
            	cnt = 0 
    if cnt!= 0:
        flen_list.append(round((cnt/fvec_len)*100, 2))
    return flen_list

           

inp_fp = None
footprint_length_dict = None
out_fp_per_orange_and_flen = None
out_fp_percentage_methylation = None
out_fp_percentage_methylation_pkl = None
# for neighborhood analysis= None
footprint_dict = None
extend_left_from_center = None
extend_right_from_center = None
strand = None
lflank = None
rflank = None

footprint_dict = pickle.load(open(sys.argv[1], "rb"))
footprint_length_dict = pickle.load(open(sys.argv[2], "rb")) 
job_id_pkl = pickle.load(open(sys.argv[3], "rb"))
 
for job_k in job_id_pkl:

    inp_fp = open(job_id_pkl[job_k][0])
    #footprint_length_dict = pickle.load(open(sys.argv[2], "rb"))
    out_fp_per_orange_and_flen = open(job_id_pkl[job_k][2], "w")
    out_fp_percentage_methylation =  open(job_id_pkl[job_k][3], "w")
    out_fp_percentage_methylation_pkl = open(job_id_pkl[job_k][4], "wb")
    # for neighborhood analysis
    #footprint_dict = pickle.load(open(sys.argv[6], "rb")) 
    extend_left_from_center = int(job_id_pkl[job_k][6]) 
    extend_right_from_center = int(job_id_pkl[job_k][7]) 
    strand = job_id_pkl[job_k][8]
    lflank = int(job_id_pkl[job_k][9])
    rflank = int(job_id_pkl[job_k][10]) 
    per_methylation_dict = defaultdict(dict)
    
    for line in inp_fp:
        # example line
        # chr3L:616610-616610^7`SRR3133326.2526820_2526820/1_overlapping`99~147	.....................	..Z..................	TACGATAATTGGTTTTTTTTT
        d_loc = [m.start() for m in re.finditer("\t", line)]
        k = line[0:d_loc[0]]
        footprint_vec = line[d_loc[0] + 1: d_loc[1]]
        per_orange = round((footprint_vec.count("F")*100)/len(footprint_vec), 3) 
        real_lengths = footprint_length_dict[k] 
        # new development: considering percentage methylation on read as feature too
        # please see page 45-47 of enhancer doc 
        m_vec = line[d_loc[1] + 1: d_loc[2]]
        
        _, left_neighbor, _ = extend_footprint_from_center(line, footprint_dict,
                                0, 0, extend_left_from_center, rflank, strand) #  
        _, right_neighbor, _ = extend_footprint_from_center(line, footprint_dict,
                                0, 0, lflank, extend_right_from_center, strand) # 
        
        per_c_l, total_upper_l, total_lower_l, total_l = get_count_and_percentage_methylation (left_neighbor) 
        per_c_r, total_upper_r, total_lower_r, total_r = get_count_and_percentage_methylation (right_neighbor) 
        if (per_c_l != "NA") and (per_c_r != "NA"):
            if float(per_c_l) <= float(per_c_r):  
                to_write_per_methylation = "\t".join ([k, str(total_upper_l), 
                                         str(total_lower_l), str(total_l), str(per_c_l)]) 
                per_methylation_dict[k]["per_c"] =  per_c_l # string data type
                per_methylation_dict[k]["total_upper"] =  total_upper_l # string data type
                per_methylation_dict[k]["total_lower"] =  total_lower_l # string data type
                per_methylation_dict[k]["total"] =  total_l # string data type
            else:
                to_write_per_methylation = "\t".join ([k, str(total_upper_r), 
                                         str(total_lower_r), str(total_r), str(per_c_r)]) 
                per_methylation_dict[k]["per_c"] =  per_c_r # string data type
                per_methylation_dict[k]["total_upper"] =  total_upper_r # string data type
                per_methylation_dict[k]["total_lower"] =  total_lower_r # string data type
                per_methylation_dict[k]["total"] =  total_r # string data type
        elif (per_c_l == "NA") and (per_c_r != "NA"): 
            to_write_per_methylation = "\t".join ([k, str(total_upper_r), 
                                     str(total_lower_r), str(total_r), str(per_c_r)]) 
            per_methylation_dict[k]["per_c"] =  per_c_r # string data type
            per_methylation_dict[k]["total_upper"] =  total_upper_r # string data type
            per_methylation_dict[k]["total_lower"] =  total_lower_r # string data type
            per_methylation_dict[k]["total"] =  total_r # string data type
        elif (per_c_l != "NA") and (per_c_r == "NA"):
            to_write_per_methylation = "\t".join ([k, str(total_upper_l), 
                                     str(total_lower_l), str(total_l), str(per_c_l)]) 
            per_methylation_dict[k]["per_c"] =  per_c_l # string data type
            per_methylation_dict[k]["total_upper"] =  total_upper_l # string data type
            per_methylation_dict[k]["total_lower"] =  total_lower_l # string data type
            per_methylation_dict[k]["total"] =  total_l # string data type
             
        else:
            to_write_per_methylation = "\t".join ([k, str(total_upper_l), 
                                     str(total_lower_l), str(total_l), str(per_c_l)]) 
            per_methylation_dict[k]["per_c"] =  "0" # string data type
            per_methylation_dict[k]["total_upper"] =  total_upper_l # string data type
            per_methylation_dict[k]["total_lower"] =  total_lower_l # string data type
            per_methylation_dict[k]["total"] =  total_l # string data type
         
        out_fp_percentage_methylation.write(to_write_per_methylation + "\n")
        percentage_orange_per_footprint = get_per_orange_for_each_footprint(footprint_vec)  
        if len(real_lengths) == 0: # Naked DNA
            to_write = k + "\t" + str(per_orange) + "\t" + str(0)  + "\t" + str(0)
            out_fp_per_orange_and_flen.write(to_write + "\n") 
        else:
      
            for l,o in zip(real_lengths, percentage_orange_per_footprint):
                to_write = k + "\t" + str(per_orange) + "\t" +  str(l) + "\t" + str(o)
                out_fp_per_orange_and_flen.write(to_write + "\n") 
    
       
    
    pickle.dump (per_methylation_dict, out_fp_percentage_methylation_pkl)
    out_fp_per_orange_and_flen.close() 
    out_fp_percentage_methylation.close()
    out_fp_percentage_methylation_pkl.close()
    inp_fp.close()
    print ("Job {} done ...".format(job_k)) 
    
