import sys
from collections import defaultdict

inp_fp = open(sys.argv[1])
out_fp = open(sys.argv[2], "w")
left_tf = {
    "T-D" : "3",
    "T-T" : "4",
    "T-N" : "5"}
right_tf = {
    "D-T" : "1",
    "T-T" : "4",
    "N-T" : "7"}
both_tf = {
    "T-T" : "4"}
accepted_labels = list (map(str, range(9)))
for f in inp_fp:
    binding_dict  = defaultdict(lambda : 0)
    tmp_f = open (f.strip())
    total_reads = 0 
    for line in tmp_f:
        binding_id = line.split("\t")[0].split("#")[-1] 
        binding_dict[binding_id] +=1
        if binding_id in accepted_labels:
            total_reads += 1 
    tmp_f.close()
    left_tf_count = 0 
    right_tf_count = 0
    both_tf_count = 0
    for k  in left_tf:
        left_tf_count += binding_dict[left_tf[k]]
    for k in right_tf: 
        right_tf_count += binding_dict[right_tf[k]]
    for k in both_tf:
        both_tf_count += binding_dict[both_tf[k]]
    
    obs_by_exp = "NA" 
    if (left_tf_count >0) and (right_tf_count >0):
        expected = (left_tf_count/total_reads)*(right_tf_count/total_reads)  
        obs_by_exp = round((both_tf_count/total_reads)/expected,3)
        to_write = "\t".join(map(str,[f.strip(), left_tf_count, right_tf_count, 
                          both_tf_count, total_reads, obs_by_exp ]))
        out_fp.write(to_write + "\n")
    else:
        to_write = "\t".join(map(str,[f.strip(), left_tf_count, right_tf_count, 
                          both_tf_count, total_reads, obs_by_exp])) 
        out_fp.write(to_write + "\n")

out_fp.close() 
inp_fp.close()
    
        
