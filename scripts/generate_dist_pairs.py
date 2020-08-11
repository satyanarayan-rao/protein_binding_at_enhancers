import os
import sys
from collections import defaultdict
from itertools import combinations
import pickle
import re
inp_fp = open(sys.argv[1])
annotation_dict = pickle.load(open(sys.argv[2],"rb"))
out_fp = open(sys.argv[3], "w")
header = "\t".join(["vplot_cl_id", "peak_id", "primary_peak_loc", 
                    "secondary_peak_loc", "chr_loc", "kmean_cl_id", "strand"])
out_bed = open(sys.argv[4], "w")
out_fp.write(header + "\n")
for f in inp_fp: 
    f_items = f.strip().split()
    v_cl_id = f_items[0]
    tmp_fp = open(f_items[1])
    for line in tmp_fp:
        # line example
        # peak_1216	0	47	53	92
        # peak_2538	0	150	62
        # <peak_id> 0 <space separated MNase peak location>  
        line_items = line.strip().split()
        peak_id = line_items[0]
        primary_peak_loc = line_items[1]
        secondary_peaks = line_items[2:]
        chr_loc = annotation_dict[peak_id]['chr_loc']
        cl_id = annotation_dict[peak_id]['cl_id'] 
        strand = annotation_dict[peak_id]['strand'] 
        for j in secondary_peaks:
            to_write = "\t".join([v_cl_id, peak_id, primary_peak_loc, j, chr_loc, cl_id, strand])
            out_fp.write(to_write + "\n")  
            chrom_details = re.split("[-:^]", chr_loc)
            chr_bed = chrom_details[0]
            chr_start = int(chrom_details[1])
            if strand == "+":
                chr_end = chr_start + int(j)
                to_write_bed = "\t".join([chr_bed, str(chr_start), str(chr_end),
                         peak_id, chr_loc, strand, primary_peak_loc, j, v_cl_id])
                out_bed.write(to_write_bed + "\n")
            elif strand == "-":
                chr_end_r = chr_start
                chr_start_r = chr_end_r - int (j) 
                to_write_bed = "\t".join([chr_bed, str(chr_start_r), str(chr_end_r), 
                         peak_id, chr_loc, strand, primary_peak_loc, j, v_cl_id])
                out_bed.write(to_write_bed + "\n")
                
    tmp_fp.close()

out_fp.close() 
inp_fp.close()
out_bed.close() 
