import os
import sys
import re
import pickle

input_fp = open(sys.argv[1])
sam_flag_dict = pickle.load(open(sys.argv[2], "rb"))
lflank = int(sys.argv[3])
rflank = int(sys.argv[4])
output_fp = open(sys.argv[5], "w")

for line in input_fp:
    # oooooooooo
    # 0123456789
    #     45
    # 
    # 15994238 - 50 - 15994009  
    delim_loc  = [m.start() for m in re.finditer("\t", line)]
    m_vec_start = int(line[delim_loc[0]+1:delim_loc[1]]) - int(line[delim_loc[6]+1:delim_loc[7]]) -  lflank
    m_vec_stop = int(line[delim_loc[0]+1:delim_loc[1]]) - int(line[delim_loc[6]+1:delim_loc[7]]) + rflank + 1
    methylation_string = line[delim_loc[10]+1:len(line)]
    read_name_str = line[delim_loc[8] + 1:delim_loc[9]]
    read_name = read_name_str[0:read_name_str.rindex("_")] 
    methylation_flank_string = methylation_string [m_vec_start:m_vec_stop]
    first_col = "`".join([line[delim_loc[2] + 1:delim_loc[3]], 
                          line[delim_loc[8] + 1:delim_loc[9]],
                          sam_flag_dict[read_name]])
    to_write= "\t".join([first_col, methylation_flank_string])
    output_fp.write(to_write + "\n")

output_fp.close() 
input_fp.close()
