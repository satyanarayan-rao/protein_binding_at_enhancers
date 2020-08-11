import os
import sys
from collections import defaultdict
import re

read_fp = open(sys.argv[1])
bed_fp = open(sys.argv[2])
out_fp = open(sys.argv[3], "w") 

read_dict = defaultdict(lambda : False)
for l in read_fp:
    read_dict[l.strip()] = True

for l in bed_fp:
   # example line: head -1 suppressed_merged_S2_to_at_50PNE_closed_peak_lf_15_rf_15_read_intersect.bed
   # chr3R	26000394	26000596	202	SRR3133326.1794773_1794773/1_overlapping`99~147	chr3R	26000595	26000625	chr3R:26000610-26000610	.	.
   d_loc = [m.start() for m in re.finditer("\t", l)]
   read_id = l[d_loc[3] + 1: d_loc[4]]
   if read_dict[read_id] == True:
       out_fp.write(l) 

out_fp.close()
bed_fp.close()
read_fp.close()
