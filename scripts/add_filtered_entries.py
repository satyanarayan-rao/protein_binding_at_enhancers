import os
import sys
import re
from collections import defaultdict

to_add_open = open(sys.argv[1]) 
existing_cnt_file = open(sys.argv[2])


# create a dictionary of exiting key: footprint length, val = count

cnt_dict = {} 
order_list = [] 
cnt = 0 
footprint_class = ""
enhancer_class = ""
for line in existing_cnt_file: 
    line_items = line.strip().split() 
    cnt_dict[line_items[0]] = int(line_items[1])
    order_list.append(line_items[0])
    if cnt == 0:
        footprint_class = line_items[3]
        enhancer_class = line_items[4]
    cnt +=1 
for l in to_add_open:
    l_items = l.strip().split()
    cnt_dict[l_items[-1]] +=1

total = 0 
for v in order_list:
    total += cnt_dict[v]
for v in order_list:
    to_print = "\t".join([v, str(cnt_dict[v]), 
                 str(round(cnt_dict[v]/total, 4)),
                 footprint_class, enhancer_class])
    print (to_print)
