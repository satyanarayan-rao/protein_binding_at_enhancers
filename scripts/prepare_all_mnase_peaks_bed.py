import sys
import re
from collections import defaultdict

enh_fp = open(sys.argv[1])
all_peak_pos_fp = open(sys.argv[2])
out_fp = open (sys.argv[3],"w")
dummy_cl_id = sys.argv[4]

mnase_peaks_in_enh = defaultdict(list)

for line in all_peak_pos_fp:
    # example line: head -1 ../vplot_lists/peak_primary_secondary.bed
    # chr2R	15518715	15518715	0	peak_5354	100	+ 
    d_loc = [m.start() for m in re.finditer("\t", line)] 
    peak_id = line[d_loc[3]+1:d_loc[4]] 
    line_s = line.strip()
    mnase_peaks_in_enh[peak_id].append(line_s)

for line in enh_fp:
    # example line: head -1 ../input_bed/reordered_20clust_147_50PNE_open_v1.bed
    # chr3R	19975479	19975979	19975740	8	-	peak_4746	38.90	SingleZ	4 

    d_loc = [m.start() for m in re.finditer("\t", line)] 
    peak_id = line[d_loc[5]+1:d_loc[6]]
    if len(mnase_peaks_in_enh[peak_id]) > 0:
        for l in mnase_peaks_in_enh[peak_id]:
            out_fp.write(l + "\n")
    else:
        if line.startswith("chr"):
            bed_line = "\t".join([
                 line[0:d_loc[0]],
                 line[d_loc[2]+1:d_loc[3]],
                 line[d_loc[2]+1:d_loc[3]],
                 "0", 
                 peak_id,
                 dummy_cl_id,
                 line[d_loc[4]+1:d_loc[5]]])
            out_fp.write(bed_line + "\n")
        else: 
            bed_line = "\t".join([
                 "chr"+line[0:d_loc[0]],
                 line[d_loc[2]+1:d_loc[3]],
                 line[d_loc[2]+1:d_loc[3]],
                 "0", 
                 peak_id,
                 dummy_cl_id,
                 line[d_loc[4]+1:d_loc[5]]])
            out_fp.write(bed_line + "\n")      
    
out_fp.close()
enh_fp.close()
all_peak_pos_fp.close() 

