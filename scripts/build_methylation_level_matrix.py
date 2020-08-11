import os
import sys
import pickle
import gzip
input_bed = open (sys.argv[1])
cpg_dict = pickle.load (open(sys.argv[2], "rb"))
gpc_dict = pickle.load (open(sys.argv[3], "rb"))
out_matrix_fp = gzip.open (sys.argv[4], "w")

for line in input_bed:
    line_items = line.strip().split()
    methylation_level_vec = []
    for loc in range(int(line_items[1]), int(line_items[2])):
        chrom_loc = "`".join([line_items[0], str(loc)])
        if (chrom_loc in cpg_dict): 
              methylation_level_vec.append (str(cpg_dict[chrom_loc]))
        elif (chrom_loc in gpc_dict):  
              methylation_level_vec.append (str(gpc_dict[chrom_loc])) 
        else:
              methylation_level_vec.append("NA")
    to_write = "\t".join(methylation_level_vec)
    to_write = "\t".join([line_items[3], to_write])
    gz_line = bytes(to_write + "\n", encoding = "ascii")
    out_matrix_fp.write (gz_line)

out_matrix_fp.close()

