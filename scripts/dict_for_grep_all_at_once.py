import os
import sys
from collections import defaultdict
import pickle 
import re
import subprocess
inp_fp = open(sys.argv[1])
out_fp = open(sys.argv[2], "wb")
chrom_to_file_dict = defaultdict(dict)
chrom_flag = defaultdict(lambda : False)
for line in inp_fp:
    # example line: head -1 select_for_site_14a_closed.sh
    # grep "chr3R:19036636-19036636\^181" flank_footprint_matrix/suppressed_merged_S2_to_all_open_and_closed_mnase_peaks_lf_15_rf_15_methylation_matrix.cluster_181.tsv | grep -w 99~147 > peak_footprints/suppressed_merged_S2_to_all_open_and_closed_mnase_peaks_cluster_181_lf_15_rf_15_methylation_matrix_for_site_peak_1341_1_sam_flag_99~147.tsv
    l_items = line.strip().split()
    chrom_and_cl = "^".join(l_items[1].split("\^"))
    if "99~147" in line:
        chrom_to_file_dict[chrom_and_cl]["99~147"] = l_items[-1]
        chrom_flag[chrom_and_cl+"99~147"] = True
    elif "83~163" in line:
        chrom_to_file_dict[chrom_and_cl]["83~163"] = l_items[-1]
        chrom_flag[chrom_and_cl+"83~163"] = True

pickle.dump(chrom_to_file_dict, out_fp)
inp_fp.close()
out_fp.close()
is_first = defaultdict(lambda : True)
sub_fp = open(sys.argv[3])
cnt_line = 0
cnt_10k = 0 
for line in sub_fp: 
    # example line: head -1 flank_footprint_matrix/suppressed_merged_S2_to_all_open_and_closed_mnase_peaks_lf_15_rf_15_methylation_matrix.cluster_181.tsv
    # chrX:317447-317447^181`SRR3133326.2448007_2448007/1_overlapping`83~163	...............................	..............z................	ATAATTTCATTTTCAAATCAATAAAAATAAA
    #d_loc = [m.start() for m in re.finditer("\t", line)] 
    chr_loc = '"' + line[0:line.find("`")] + '"'
    if "99~147" in line: 
        if chrom_flag[chr_loc+"99~147"] == True:
            fname = chrom_to_file_dict[chr_loc]["99~147"]
            l_s = line.strip() 
            l_s = l_s.replace("`", "\\`").replace("\t", "\\\t")
            if (os.path.isfile(fname) == True) and (os.path.getsize(fname) !=0) and (is_first[chr_loc+"99~147"] == True):
                cmd = " ".join(["echo", l_s, ">", fname])
                os.system(cmd)
                is_first[chr_loc+"99~147"] = False
            else:
                cmd = " ".join(["echo", l_s, ">>", fname])
                os.system(cmd)
                is_first[chr_loc+"99~147"] = False
                   
            
    elif "83~163" in line:
        if chrom_flag[chr_loc+"83~163"] == True:
            fname = chrom_to_file_dict[chr_loc]["83~163"]
            l_s = line.strip() 
            #l_s = l_s.replace("`", "\\`")
            l_s = l_s.replace("`", "\\`").replace("\t", "\\\t")
            if (os.path.isfile(fname) == True) and (os.path.getsize(fname) !=0) and (is_first[chr_loc +"83~163"] == True):
                cmd = " ".join(["echo", l_s, ">", fname])
                os.system(cmd)
                is_first[chr_loc+"83~163"] = False
                
            else:
                cmd = " ".join(["echo", l_s, ">>", fname])
                os.system(cmd)
                is_first[chr_loc+"83~163"] = False
    cnt_line +=1 
    if cnt_line % 10000 == 0:
        cnt_10k += 1
        print ("completed {} lines ...".format(cnt_10k*10000))
        
 
sub_fp.close()
