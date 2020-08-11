import os
import sys
from collections import defaultdict, OrderedDict
import pandas as pd
import sys
import gzip
import getopt
import re
import math
import numpy as np
import scipy
from scipy.signal import find_peaks
from numpy import pi
from scipy.signal import savgol_filter
import json

def count_footprints (vec_list):
    """
    vec_list ~ [ ".......FFFFF.....FF....FFFF", "FFFFFFFFFF....FFFFF"] 
    """
    total_bases = 0 
    total_footprints = 0 
    for vec  in vec_list: 
        total_footprints += vec.count("F")
        total_bases += len(vec)
    fraction = round(total_footprints/total_bases, 3)
    return total_bases, total_footprints, fraction
# sys.argv[1]: take clustered footprint file
# sys.argv[2]: take kde file
footprint_fp = open(sys.argv[1])
#kde_fp = open(sys.argv[2])
out_fp = open(sys.argv[3], "w")
out_txt_fp = open(sys.argv[4], "w")
percentage_molecule_fp = open(sys.argv[5], "w")
peak_label  = sys.argv[6]

all_cls = {}
footprint_dict = defaultdict(list) 
total_dna_molecules = 0
for line in footprint_fp: 
    line_items = line.strip().split()
    cl_id = int (line_items[0].split("#")[-1])
    all_cls[cl_id] = True
    footprint_dict[cl_id].append(line_items[1])
    total_dna_molecules += 1

cl_ids_sorted = sorted(list(all_cls.keys()))
reassingned_cl_footprint_dict = defaultdict(list)   

for i in range (1, len(cl_ids_sorted) + 1): 
    reassingned_cl_footprint_dict[i] = footprint_dict[cl_ids_sorted[i-1]] 

stats_dict = defaultdict(dict)
naked_dna_cluster = None
min_fraction = 100 
for k in reassingned_cl_footprint_dict: 
    total_bases, total_footprints, fraction  =  count_footprints(reassingned_cl_footprint_dict[k])
    if min_fraction > fraction: 
        min_fraction = fraction
        naked_dna_cluster = k

    stats_dict[k]["total"] = total_bases
    stats_dict[k]["footprints"] = total_footprints
    stats_dict[k]["fraction"] = fraction
    stats_dict[k]["total_molecules"] = len(reassingned_cl_footprint_dict[k])  
    
#nuc_range_start = 120
#kde_sum_dict = defaultdict(lambda : 0 )

kde_df = pd.read_csv(sys.argv[2], header = None, sep = "\s+")

kde_df.columns = ["x_tics", "kde", "cl_id"]
# I have to find peaks with different values of mean and standard deviation means
peak_dict =  defaultdict(lambda : defaultdict(OrderedDict))
#mean_ratio_to_try = [0.4, 0.5, 0.6, 0.75, 1] 
mean_ratio_to_try = [0] 
all_clusters = sorted(kde_df.cl_id.unique())

#header_str = ["cl_id", "ratio", "x_tics", "peak_loc", ""]
#header_str = "%6s\t%10s\t%15s\t%20s\t%15s\t%10s\t%10s\t%15s\t%15s\t%15s\t%10s"%("cl_id", "ratio", 
#                    "peak_loc", "peak_y", "total_peaks", "mean_val", "height",
#                    "per_orange", "footprints", "total", "prob_tf") 
header_str = "%15s\t%6s\t%15s\t%10s\t%20s\t%20s\t%10s\t%15s\t%15s\t%15s\t%10s"%(
            "enhancer_id", "cl_id", "per_orange", "gt_one_peak", "peak_pos", 
            "peak_y", "total_peaks", "mean_kde", "footprints", "total", "prob_tf")
out_txt_fp.write(header_str + "\n")
nuc_cluster = None
prob = 2 
for cl in all_clusters: 
    kde_vals = np.array(kde_df.loc[kde_df.cl_id == cl, ]["kde"]) 
    mean_val = kde_vals.mean() 
    ####### calculate area under the density curve ##### 
    kde_x_tics = np.array(kde_df.loc[kde_df.cl_id == cl, ]["x_tics"]) 
    total_area = np.trapz(kde_vals, x = kde_x_tics) # No need to give dx - because they are evenly spaced
    area_le_100 = np.trapz(kde_vals[kde_x_tics<=100], x  = kde_x_tics[kde_x_tics<=100])
    prob_tf = round(area_le_100/total_area, 2) 
    if prob > prob_tf: 
        prob = prob_tf 
        nuc_cluster = cl
    for r in mean_ratio_to_try:
        peaks_x_tics = [round (i,2) for i in kde_df.loc[find_peaks(kde_vals)[0],]["x_tics"].tolist() ] 
        peaks_y_vals = [str(format(i, "0.2E")) for i in kde_df.loc[find_peaks(kde_vals)[0],]["kde"].tolist() ] 
        str_peak = "`".join(map(str, peaks_x_tics))
        y_vals_at_peak = "`".join(peaks_y_vals)
        gt_one_peak = 0
        if len(peaks_x_tics) > 1: 
            gt_one_peak = 1
        peak_dict[str(cl)][str(r)]["peak_loc"] = "`".join(map (str, peaks_x_tics))
        peak_dict[str(cl)][str(r)]["peak_y"] = y_vals_at_peak 
        peak_dict[str(cl)][str(r)]["total_peaks"] = len(peaks_x_tics)
        peak_dict[str(cl)][str(r)]["mean_val"] =  mean_val
        #peak_dict[str(cl)][str(r)]["height"] =  r*mean_val  # no longer used 
        peak_dict[str(cl)][str(r)]["percent_orange"] = stats_dict[cl]["fraction"]
        peak_dict[str(cl)][str(r)]["total"] = stats_dict[cl]["total"]
        peak_dict[str(cl)][str(r)]["footprints"] = stats_dict[cl]["footprints"]
        to_write = "%15s\t%6s\t%15s\t%10s\t%20s\t%20s\t%10s\t%15s\t%15s\t%15s\t%10s"%( 
                    peak_label, 
                    cl, 
                    round(stats_dict[cl]["fraction"]*100, 2), 
                    gt_one_peak,
                    str_peak, 
                    y_vals_at_peak,   
                    len(peaks_x_tics), 
                    format(mean_val, "0.2E"), 
                    #str(format(r*mean_val, "0.2E")),
                    stats_dict[cl]["footprints"], 
                    stats_dict[cl]["total"], 
                    str(prob_tf)) 
        out_txt_fp.write(to_write + "\n")
        

# Now we have to decide cluster associactions: naked DNA, TF and NUC 
# since we have already decided naked DNA cluster based on %orange 
# see `naked_dna_cluster` variable
naked_and_nuc_cluster = [naked_dna_cluster, nuc_cluster]
tf_cluster = None 
for k in all_clusters:
    if k in naked_and_nuc_cluster: 
        continue
    else: 
        tf_cluster = k 

print (all_clusters)
print ([naked_dna_cluster, nuc_cluster, tf_cluster])
# calculate percentage DNA molecules in each cluster 
percentage_naked = round(stats_dict[naked_dna_cluster]["total_molecules"]/total_dna_molecules, 2) * 100
percentage_tf = round(stats_dict[tf_cluster]["total_molecules"]/total_dna_molecules, 2) * 100
percentage_nuc = round(stats_dict[nuc_cluster]["total_molecules"]/total_dna_molecules, 2) * 100

header = "%10s\t%20s\t%10s\t%20s\t%10s\t%20s\t%10s"%(
          "naked_cl", "percentage_naked", "tf_cl", "percentage_tf", "nuc_cl",
          "percentage_nuc", "total_molecules")
percentage_molecule_fp.write(header + "\n")
vals = "%10s\t%20s\t%10s\t%20s\t%10s\t%20s\t%10s"%(
         naked_dna_cluster, percentage_naked, tf_cluster, percentage_tf, 
         nuc_cluster, percentage_nuc, total_dna_molecules)
percentage_molecule_fp.write(vals+"\n")
tt = dict(peak_dict)
json.dump (tt, out_fp) 
out_fp.close()
out_txt_fp.close()
percentage_molecule_fp.close()
