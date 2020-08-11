import os
import sys
label_fp = open(sys.argv[1])
footprint_fp = open(sys.argv[2])
out_fp = open(sys.argv[3], "w")
label_dict = {} 
for line in label_fp:
    line_items = line.strip().split() 
    label_dict[line_items[0]] = line_items[1] 
for line in footprint_fp:
    k = line[0:line.find("\t")]
    to_write = "\t".join([k + "#" + label_dict[k], line[line.find("\t") + 1:]])
    out_fp.write(to_write) 

label_fp.close()
footprint_fp.close()
out_fp.close()
