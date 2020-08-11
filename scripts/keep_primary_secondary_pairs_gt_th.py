import sys
from collections import defaultdict

inp_fp = open(sys.argv[1])
min_dist = int(sys.argv[2])

peak_dict = defaultdict(lambda : False)
peak_list = [] 
peak_dist = defaultdict (list)
for line in inp_fp:
    l_items = line.strip().split() 
    peak_id = l_items[0] 
    cl_id = l_items [1] 
    dist = int(l_items[2])
    if dist < min_dist:
        continue
    else:
        if peak_dict[peak_id] == False:
            peak_list.append(peak_id)
            peak_dict[peak_id] = True
            peak_dist[peak_id].append(dist)
            print ("\t".join([peak_id, cl_id, str(dist)]))
            continue
        
        if len(peak_dist[peak_id]) >=1:
            is_lt = False 
            for d in peak_dist[peak_id]:
                if dist - d < min_dist:
                    is_lt = True
            if is_lt == False:
                peak_dist[peak_id].append(dist) 
                print ("\t".join([peak_id, cl_id, str(dist)]))

#for k in peak_list:
#    for d in 
                
