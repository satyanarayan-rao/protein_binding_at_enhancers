import sys 
import pandas as pd

dt = pd.read_csv(sys.argv[1], sep = "\s+", header = "infer") 

#file_name_without_sam_flag	total_reads	total_naked	total_tf	total_nuc	occup_naked	occup_tf	occup_nuc	total_discard	occup_discard
#binding_labeled_reads_for_peak/suppressed_merged_S2_to_all_open_and_closed_mnase_peaks_cluster_181_lf_15_rf_15_site_peak_2330_1_binding_label.cond8.tsv	37	2	1	32	5.41	2.7	86.49	2	5.41
print(dt.columns)
cols_list = list(dt.columns)
renamed_list = [] 
for c in cols_list:
    if c in ['occup_naked', 'occup_tf', 'occup_nuc', 'total_reads']:
        r = "with_discard_" + c
        renamed_list.append(r)
    else:
        renamed_list.append(c)

dt.columns = renamed_list 

dt["total_reads"] = dt.with_discard_total_reads - dt.total_discard
dt["occup_naked"] = (dt.total_naked/dt.total_reads)*100
dt["occup_tf"] = (dt.total_tf/dt.total_reads)*100
dt["occup_nuc"] = (dt.total_nuc/dt.total_reads)*100

dt.to_csv(sys.argv[1] + ".caliberated.tsv", float_format = "%0.2f", sep = "\t", index = False)
