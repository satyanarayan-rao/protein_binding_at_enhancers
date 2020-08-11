import os
import sys
import pickle
read_samflag_dict = {}
out_fp = open(sys.argv[1], "wb")
for l in sys.stdin:
    read1 = l 
    read2 = next(sys.stdin)
    tab_loc = read1.find("\t")
    s_flag_r1 = read1[tab_loc + 1: len(read1) - 1] 
    s_flag_r2 = read2[tab_loc + 1: len(read2) - 1]
    read_samflag_dict[read1[0:tab_loc]] = s_flag_r1 + "~" + s_flag_r2
pickle.dump(read_samflag_dict, out_fp)
out_fp.close()
