import os
import sys
import re
from collections import defaultdict
from collections import OrderedDict
import math
# sys.argv[1] = bed file with last field methylated reads
# sys.argv[2] = output bed file with defined footprints

def call_footprints(methylation_vec): 
    # Example: 
    # -------------------------
    
    # % : unmethylated  C 
    # $ : methylated C 
    # . : non-C
    # ^ : footprint
    # Q : Not allowed  
    # Rule for calling fooptrint: if unmethylated cytosine is found 
    # between two methylated cytosines (+- 2bp buffer; owing to the fact that 
    # methyltransferases methylate in CG, CHG range)  
    
    # Original read:  $ . $ . $ . . . $ . . . . % . . . . . $ . . . % . . . $ . . . .
    #                 0 1 2 3 4 5 6 7 8 9 10111213141516171819202122232425262728293031
    # not allowed:    Q Q Q Q Q Q Q Q Q Q Q . . % . . . Q Q Q Q Q . % . Q Q Q Q Q Q Q 
    #                 . . . . . . . . . . | . . % . . . | . . . | . % . | . . . . . .
    # Footprints:     . . . . . . . . . . . F F F F F F . . . . . F F F . . . . . . .
    
    cnt = 0
    footprint_vec = []
    footprint_chr_list = []
    upper_case_loc = []
    lower_case_loc = []
    footprint_vec = []
    to_return_vec = ""
    footprint_length_dict = OrderedDict()
    vec_len = len(methylation_vec)
    
    for i in range(len(methylation_vec)):  # find all the positions of upper case
        l = methylation_vec[i]
        if (l == l.upper()) and (l != ".") : 
            upper_case_loc.append(i)
        elif (l == l.lower()) and (l != "."): 
            lower_case_loc.append(i) 

    not_allowed_pos_vec = list(methylation_vec)
    
    if len (upper_case_loc)>=2: # check if at least two upper case letters are there in the string
        for l in upper_case_loc:
            start_x = max(0, l - 2) # 4 - 2 = 2
            end_x = min (l+3, vec_len) # 4 + 3 = 7 # will fill x from 2 to 6 for `$` at 4
            for x_pos in range(start_x, end_x): 
                not_allowed_pos_vec[x_pos] = "Q" 
        footprint_reduced_vec = not_allowed_pos_vec
        boundary_locations = [] 
        #print("".join(footprint_reduced_vec))    
        #print("\n#####\n")
        next_set = False 
        for i in range(vec_len - 1):
            if next_set == True: 
                next_set = False
                continue
            li = not_allowed_pos_vec[i]
            li_plus_one = not_allowed_pos_vec[i+1]
            if (li == li.upper()) and (li != ".") and (li_plus_one in [".", li_plus_one.lower()]): # handle `Q.`;
                footprint_reduced_vec[i] = "|"
                boundary_locations.append(i)
            elif (li in [".", li.lower()]) \
                 and (li_plus_one == li_plus_one.upper()) \
                 and (li_plus_one != "."): # handle .Q
                 footprint_reduced_vec[i] = li
                 footprint_reduced_vec[i+1] = "|"
                 next_set = True
                 boundary_locations.append(i+1)
            elif (li == ".") or (li == li.lower()):
                footprint_reduced_vec[i] = li
            else: 
                footprint_reduced_vec[i] = "."
        last_letter = footprint_reduced_vec[vec_len - 1] 
        if (last_letter == last_letter.upper()) and (last_letter != "."):
             footprint_reduced_vec[vec_len - 1] = "." 
         
        #print ("".join(footprint_reduced_vec))    
          
        footprint_chr_list.extend("."*boundary_locations[0]) 
        for pipe_loc in range(len(boundary_locations) - 1):
            current_loc = boundary_locations[pipe_loc]
            next_loc = boundary_locations[pipe_loc + 1 ]
            unmeth_str = "".join(footprint_reduced_vec[current_loc + 1 : next_loc])
            unmeth_loc = []
            for l in unmeth_str: 
                if l == ".":
                    continue
                elif l == l.lower():
                    unmeth_loc.append(l)
            if len(unmeth_loc) >=1:
                footprint_chr_list.extend(".")
                footprint_chr_list.extend("F"*(next_loc - current_loc - 1 ))
                footprint_length_dict["-".join([str(current_loc), str (next_loc)])] =  next_loc - current_loc - 1
            else: 
                footprint_chr_list.extend(".")
                footprint_chr_list.extend("."*(next_loc - current_loc - 1 ))
        footprint_chr_list.extend("."*(vec_len - boundary_locations[-1])) 
                
        footprint_vec = "".join(footprint_chr_list) 
        return footprint_vec, footprint_length_dict 
        
    else:
       to_return_vec = "."*len(methylation_vec) 
       return to_return_vec, footprint_length_dict

input_fp = open(sys.argv[1])
out_fp = open(sys.argv[2], "w") 
footprint_length_fp = open(sys.argv[3], "w")

for line in input_fp:
    d_loc = [m.start() for m in re.finditer("\t", line)]
    methylation_vec = line[d_loc[-1] + 1: len(line) - 1]
    footprint_vec, footrprint_length_dict = call_footprints(methylation_vec) 
    to_write = line[0:d_loc[-1]] + "\t" + footprint_vec
    out_fp.write(to_write + "\n")
     
    for k in footrprint_length_dict:
        to_write = "\t".join([ line[0:d_loc[-1]], k, str(footrprint_length_dict[k]) ])
        footprint_length_fp.write(to_write + "\n")
out_fp.close()
footprint_length_fp.close()
