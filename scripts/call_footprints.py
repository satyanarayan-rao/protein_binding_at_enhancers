import os
import sys
import re
from collections import defaultdict
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
    
    # Original read:  $ . $ . $ . . . $ . . . . % . . . . . $ . . . % . . . $ . . . .
    #                 0 1 2 3 4 5 6 7 8 9 10111213141516171819202122232425262728293031
    # Footprint read: . . . . . . . . . . ^ ^ ^ ^ ^ ^ ^ . . . . ^ ^ ^ ^ ^ . . . . . . 
    #                 0 1 2 3 4 5 6 7 8 9 10111213141516171819202122232425262728293031
    cnt = 0
    footprint_vec = []
    upper_case_loc = []
    footprint_vec = []
    to_return_vec = ""
    footrprint_length_dict = {}
    for i in range(len(methylation_vec)): 
        l = methylation_vec[i]
        if (l == l.upper()) and (l != ".") : 
            upper_case_loc.append(i)


    if len (upper_case_loc)>=2:
        footprint_vec.extend("."*upper_case_loc[0]) # first append `.` untill the first capital letter
        for i in range(len(upper_case_loc) - 1):
            current_loc = upper_case_loc[i]
            next_loc = upper_case_loc[i + 1]
            string_in_between = methylation_vec[current_loc + 1: next_loc] # find footrpint in this
            unmethylated_loc = []
            for j in range(len(string_in_between)): 
                l = string_in_between[j]
                if (l == l.lower()) and (l != "."): 
                    unmethylated_loc.append(current_loc + j)
                else:
                    continue
            if len (unmethylated_loc) >=1: 
                left_boundary = (unmethylated_loc[0] + current_loc)//2
                right_boundary = (unmethylated_loc[-1] + next_loc)//2
                footprint_vec.append("."*(left_boundary - current_loc)) 
                footprint_vec.append ("F"*(right_boundary - left_boundary + 1 )) 
                footprint_vec.append("."*(next_loc - right_boundary - 1))
                footrprint_length_dict["-".join([str(left_boundary), str(right_boundary)])] = right_boundary - left_boundary + 1
            else:
                footprint_vec.append("."*(next_loc - current_loc))
        # append the last one
        footprint_vec.append("."*( len(methylation_vec) - upper_case_loc[-1]))
        to_return_vec = "".join(footprint_vec)
      
        return to_return_vec, footrprint_length_dict 
                
    else:
       to_return_vec = "."*len(methylation_vec) 
       return to_return_vec, footrprint_length_dict

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
