import sys
from collections import defaultdict
import re
import pickle
from itertools import combinations

inp_fp = open (sys.argv[1])
th = int(sys.argv[2])
chr_loc_dict = defaultdict(list)
peak_list_dict = defaultdict(list)
peak_chr = {} 
peak_list = [] 
is_peak_in_list = defaultdict(lambda : False)
for line in inp_fp:
    # example line: head -1 all_open_mnase_peaks_peak_id_chr_loc_cl_id.tsv 
    # peak_4746_1	chr3R:19975666-19975666^199	199	+
    l_items = line.strip().split() 
    peak_id_without_offset = "_".join(l_items[0].split("_")[0:2]) 
    peak_list_dict[peak_id_without_offset].append(l_items[0])
    chr_loc = int(l_items[1].split("-")[0].split(":")[1])
    chr_loc_dict[peak_id_without_offset].append(chr_loc)
    peak_chr[peak_id_without_offset] = l_items[1].split("-")[0].split(":")[0]
    if is_peak_in_list[peak_id_without_offset] == False:
        peak_list.append(peak_id_without_offset)
        is_peak_in_list[peak_id_without_offset] = True
for peak in peak_list:
    tot_peaks = len(peak_list_dict[peak])
    if tot_peaks > 1: 
        for i, j in combinations(range(tot_peaks), 2):
            if abs(chr_loc_dict[peak][i] - chr_loc_dict[peak][j]) >= th:
                if chr_loc_dict[peak][j] > chr_loc_dict[peak][i]: 
                    print ("\t".join([peak_chr[peak],
                                      str(chr_loc_dict[peak][i]), 
                                      str(chr_loc_dict[peak][j]),
                                      peak_list_dict[peak][i] + "@" + peak_list_dict[peak][j],
                                      peak_chr[peak] + ":" + str(chr_loc_dict[peak][i]) + "-" + str(chr_loc_dict[peak][i]) + "^199", 
                                      "+",
                                      "0", 
                                     str(chr_loc_dict[peak][j] - chr_loc_dict[peak][i]), 
                                     "199"]))
                else:
                    print ("\t".join([peak_chr[peak],
                                      str(chr_loc_dict[peak][j]), 
                                      str(chr_loc_dict[peak][i]),
                                      peak_list_dict[peak][j] + "@" + peak_list_dict[peak][i],
                                      peak_chr[peak] + ":" + str(chr_loc_dict[peak][j]) + "-" + str(chr_loc_dict[peak][j]) + "^199", 
                                      "+",
                                      "0", 
                                     str(chr_loc_dict[peak][i] - chr_loc_dict[peak][j]), 
                                     "199"]))
    else:
        continue
inp_fp.close() 
