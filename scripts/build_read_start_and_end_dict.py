import sys
import re
from collections import defaultdict
import pickle
inp_fp = open(sys.argv[1])
out_fp = open(sys.argv[2], "wb")

read_dict = defaultdict(dict)
for line in inp_fp:
    # example line: head -1 footprints/suppressed_merged_S2_to_all_peaks_lf_15_rf_15_tf_footprint.bed
    # chr3R	19975740	19975741	chr3R:19975740-19975740^8	.	-	chr3R	19975665	19975942	SRR3133326.2817484_2817484/1_overlapping`83~163	.	......................................................................................................................................................................................................................................................................................	...........z...........H.....X.....Z......Z.Z.............H................Z.....H.Z.....H..........H.Z..............H....Z.............H.......H...........z.............z...z.....z.........z.............h.....z......h..........z......zx.................h..................z....	AATAACATTTCAATAATTCATATGCTTCAGCAAACGAATCTCGCGATATATACCCTTAGCATCCTCCTCCCTCTCGAACGAGCGCATAAGCCTTTTCATAGCGAAATAATTTTCAATGCCTCGTAATCTCACCTTAGCCACCTAGCCAAAAAAACCACCACCCAAAAAACAAACAAACTCATATATATCCAAAAACTCCCAAAAACTTTCATTAATTACCACTCTCACAAACTCCAACATTCTAAATTAAATATACAAAAATATATCTTCTACAATAA
    d_loc = [m.start() for m in re.finditer("\t", line)] 
    chrom = line[d_loc[5] + 1: d_loc[6]]
    read_start = int(line[d_loc[6] + 1: d_loc[7]])
    read_end = int(line[d_loc[7] + 1: d_loc[8]])
    read_id = line[d_loc[8] + 1: d_loc[9]] 
    read_dict[read_id]["chrom"] = chrom
    read_dict[read_id]["start"] = read_start
    read_dict[read_id]["end"] = read_end
    
pickle.dump(read_dict, out_fp)
out_fp.close() 
inp_fp.close()
