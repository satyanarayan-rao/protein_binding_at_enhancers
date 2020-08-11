import os
import sys
import pickle
import re
bed_fp = open(sys.argv[1])
flag_dict = pickle.load(open(sys.argv[2], "rb"))
out_fp = open(sys.argv[3], "w") 


for line in bed_fp:
    d_loc = [m.start() for m in re.finditer("\t", line)] 
    read_name_str = line[d_loc[8] +1 : d_loc[9]]
    read_name = read_name_str[0:read_name_str.rindex("_")] 
    sam_flag = flag_dict[read_name] 
    new_line = line[0:d_loc[9]] +  "`" + sam_flag + line[d_loc[9]:]
    out_fp.write(new_line)
out_fp.close()
bed_fp.close()
