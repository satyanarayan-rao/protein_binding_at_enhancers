import os
import sys
import re
from collections import defaultdict
inp_fp = open(sys.argv[1]) 
out_fp = open(sys.argv[2], "w")


def capital_perentage_and_stretch_of_unbound(fvec, mvec):
    len_vec = len(fvec)
    check_vec = "."*len_vec
    total_capital = 0 
    total_small = 0
    left_extreme = 0 
    right_extreme = len_vec 
    first_found = False
    base_counter = -1
    last_capital = -1   
    for c in mvec: 
        base_counter +=1
        if c == ".":
            continue 
        elif c == c.lower(): 
            total_small +=1 
        elif c == c.upper():
             total_capital +=1 
             if first_found == False: 
                 left_extreme = base_counter 
                 first_found = True
             else:
                 last_capital = base_counter
        else:
            continue 
    naked_dna_stretches = [] 
    if first_found == True: # there was at least one methylation event
         if last_capital == -1: # there was only one capital
             naked_dna_stretches.append(left_extreme) 
             naked_dna_stretches.append(0)
             naked_dna_stretches.append (right_extreme - left_extreme)
         else:
             naked_dna_stretches.append(left_extreme) 
             naked_dna_stretches.append(last_capital - left_extreme - 1)
             naked_dna_stretches.append(right_extreme - last_capital - 1)  
    else:
        naked_dna_stretches = [0, 0, 0]
    total = total_capital + total_small 
    percentage_capital = "NA"
    if total > 0:
    	percentage_capital = round(total_capital/total, 2)
    
    return len_vec, total, percentage_capital, naked_dna_stretches 
         
        

for line in inp_fp:
    # chr3R	26000913	26000914	chr3R:26000913-26000913	.	.	chr3R	26000649	26000928	SRR3133326.2348423_2348423/1_overlapping`83~163	.	........................................................................................................................................................................................................................................................................................	.....Z..H.Z....Z..............H.ZX...............................z..........................h....x...........z..h.z.......x......h........z..............h.z........................h.zx.z....x...........h..........h......................z...........z.......h...z.....z.............	ATTTCGAAGCGATTCGAATTAGTTCGTATAGCGGCTCTCTATTATAAATACATTTCCCCATCCTCAAATTCTATTTTATTTTCTTTTATTTTACTCTACTTCACTACTCAATACATCTCTCTACCAAAAACCACCTCCAAAAAATATCCTCTTACACTAAAATCACCTACATTTTCCTTAACAACAAACAACATACATATATACCAACTTATAACATATATATATATTCTATTCACACCTATATACTCAATATTTTACTCACCTACATTTCTTTTTCTGT
    d_loc = [m.start() for m in re.finditer("\t", line)]
    fvec = line[d_loc[10] + 1: d_loc[11]] 
    mvec = line[d_loc[11] + 1: d_loc[12]]
    read_id = line[d_loc[8] + 1: d_loc[9]] 
    peak_id = line[d_loc[2] + 1: d_loc[3]] 
    check_vec = "."*len(fvec)
    if check_vec == fvec:
        len_vec, total, pcap, naked_dna_stretches = capital_perentage_and_stretch_of_unbound(fvec, mvec) 
        to_write = "\t".join([
                    read_id, peak_id, str (len_vec),
                    str(total), str(pcap), 
                    "\t".join(map(str,naked_dna_stretches)),
                    str(max(naked_dna_stretches)),
                    str(max(naked_dna_stretches[0], naked_dna_stretches[2]))])
        out_fp.write(to_write + "\n")
    else:
        continue

out_fp.close()
inp_fp.close()
