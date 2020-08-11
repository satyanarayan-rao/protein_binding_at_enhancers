import os
import sys
from collections import defaultdict
import re
import pandas as pd
def get_start_and_end(vec):
    len_vec = len(vec)
    st_cnt = 0 
    en_cnt = 0
    for i in range(len_vec):
        if vec[i] == "M":
            st_cnt +=1 
        else:
            break
    for i in range(len_vec - 1, -1, -1):
        if vec[i] == "M":
            en_cnt +=1 
        else:
            break
    return st_cnt, en_cnt 
def add_edge_colors(vec, st, en):
    vec_l = list(vec)
    for i in range (st, en + 1):
        vec_l[i] = 'E'
    ret_vec = "".join(vec_l)
    return ret_vec




inp_fp = open(sys.argv[1])
out_file = sys.argv[2]

start = []
end = []
read_id_list = []
read_label_list = [] 
fp_list = [] 
for line in inp_fp:
    hash_loc = line.find("#")
    tab_loc = line.find("\t")
    #print ([hash_loc, tab_loc])
    read_id = line[0:hash_loc]
    label = line[hash_loc + 1: tab_loc] 
    fp_str = line[tab_loc + 1: len(line) - 1]
    fp_start, fp_end = get_start_and_end(fp_str)  
    start.append(fp_start)
    end.append(fp_end)
    read_id_list.append(read_id)
    read_label_list.append(label) 
    fp_list.append(fp_str)

df_dict = {"read_id": read_id_list, "read_label" : read_label_list, 
                  "st": start, "en" : end, "fp_str" : fp_list }
df = pd.DataFrame (df_dict)
#print(df.head())

sorted_df = df.sort_values(by = ["read_label", "st", "en"], ascending = [True, True, False])  

sorted_df["tagged"] = sorted_df["read_id"] + "#" + sorted_df["read_label"] 

sorted_df[["tagged", "fp_str"]].to_csv(out_file, sep = "\t", header = False, index = False)

