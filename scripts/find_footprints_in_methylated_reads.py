import os
import sys
from collections import defaultdict
import re
import math
# sys.argv[1] = bed file with last field methylated reads
# sys.argv[2] = output bed file with defined footprints
# sys.argv[3] = threshold to define footprint

def find_fotprints_on_vec (methylation_vec, th):
    # Example: 
    # -------------------------
    
    # % : unmethylated  C 
    # $ : methylated C 
    # . : non-C 
    # threshold = 3 
    
    # Original read:  $.$.$...$....%.....$...%...$....
    # Footprint read: ........$....%.....$...%...$....
    # Footprints:     --------+---f1-----+--f2---+----
    
    # threshold = 4: second one will not be called
    # Original read:  $.$.$...$....%.....$...%...$....
    # Footprint read: ........$....%.....$............
    # Footprints:     --------+---f1-----+------------
    # in actual read: 
    #         small letter: unmethylated
    #       capital letter:   methylated  
    #                  dot:        non-C
    cnt = 0
    footprint_vec = []
    upper_case_loc = []
    footprint_vec = []
    to_return_vec = ""
    for i in range(len(methylation_vec)): 
        l = methylation_vec[i]
        if (l == l.upper()) and (l != ".") : 
            upper_case_loc.append(i)
    # from the beggining of methylation_vec to the very first Capital letter
    if len (upper_case_loc)>=2:
        footprint_vec.extend("."*upper_case_loc[0]) 
        previous_was_footprint = False
        for i in range(len(upper_case_loc) - 1):
            current_loc = upper_case_loc[i]
            next_loc = upper_case_loc[i + 1]
            string_in_between = methylation_vec[current_loc + 1: next_loc] # find footrpint in this
            unmethylated_loc = [] 
            if next_loc - current_loc - 1 >= th: 
                #for j in range(len(string_in_between)): 
                #    l = string_in_between[j]
                #    if (l == l.lower()) and (l != "."): 
                #        unmethylated_loc.append(current_loc + j)
                #    else:
                #        continue
                #if len(unmethylated_loc) == 0:
                #    # footprint not found between i and i+1; only `.` was found in between
                #    # apending in footprint_vec 
                #    footprint_vec.extend(methylation_vec[current_loc: next_loc])
                #    footprint_vec[current_loc] = "." # have to reduce to non-footprint
                #else:
                #    leftmost_unmethylated = unmethylated_loc[0]
                #    rightmost_unmethylated = unmethylated_loc[-1]
                #    print(next_loc, leftmost_unmethylated, rightmost_unmethylated, next_loc - rightmost_unmethylated - 1)
                #    if (leftmost_unmethylated - current_loc - 1 >= th) and \
                #       (next_loc - rightmost_unmethylated - 1 >= th): 
                #        #  found footprint 
                #        footprint_vec.extend(methylation_vec[current_loc: next_loc])
                #    else:
                #        footprint_vec.extend("."*(next_loc - current_loc))
                footprint_vec.extend (methylation_vec[current_loc: next_loc])
                previous_was_footprint = True
            elif previous_was_footprint == True:
                footprint_vec.append(methylation_vec[current_loc])
                footprint_vec.extend("."*(next_loc - current_loc - 1))
                previous_was_footprint = False
            else:
                footprint_vec.extend("."*(next_loc - current_loc))
        if footprint_vec[upper_case_loc[-2]] != methylation_vec[upper_case_loc[-2]]:
            footprint_vec.extend(methylation_vec[upper_case_loc[-1]: len(methylation_vec)]) 
            footprint_vec[upper_case_loc[-1]] = "." 
        else: 
            footprint_vec.extend(methylation_vec[upper_case_loc[-1]: len(methylation_vec)]) 
        to_return_vec = "".join(footprint_vec)
        
    else: 
       to_return_vec = "."*len(methylation_vec) 
       return to_return_vec
    return to_return_vec
input_fp = open(sys.argv[1])
out_fp = open(sys.argv[2], "w")
th = int(sys.argv[3])
for line in input_fp:
    d_loc = [m.start() for m in re.finditer("\t", line)]
    methylation_vec = line[d_loc[-1] + 1: len(line) - 1]
    footprint_vec = find_fotprints_on_vec(methylation_vec,
                       th = th) 
    to_write = line[0:d_loc[-1]] + "\t" + footprint_vec
    out_fp.write(to_write + "\n")

out_fp.close()
