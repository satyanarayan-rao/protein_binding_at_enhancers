import sys
import re
from collections import defaultdict
import pickle
import numpy as np
import random
import math
from scipy import stats

inp_fp_83 = open(sys.argv[1])
inp_fp_99 = open(sys.argv[2])
out_fp = open(sys.argv[4], "w")
out_fp_verb = open(sys.argv[4] + "_verbose.tsv", "w")

observed_label_list = [] 
observed_label_count_dict = defaultdict(lambda : 0)
for line in inp_fp_83:
    # example line: head -1 /beevol/home/satyanarr/workplace/projects/enhancer-cooperativity_smk/annotated_cobinding_status/suppressed_merged_S2_to_all_open_and_closed_mnase_peaks_cluster_199_lf_15_rf_15_sec_peak_at_52_lex_300_rex_300_site_peak_1027_4_sam_flag_83~163_cobinding.tsv
    # chr3L:3385499-3385499^199`SRR3133326.23743952_23743952/1_overlapping`83~163#0	MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMFF.FFFFFFFFFFFFFFFFFFFFFF............FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF...........................................................................................................FFFFFFFFFFFFFFFFMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
    label = line[0:line.find("\t")].split("#")[-1]
    observed_label_list.append(label)
    observed_label_count_dict[label] +=1 
for line in inp_fp_99:
    # example line:  head -1 /beevol/home/satyanarr/workplace/projects/enhancer-cooperativity_smk/annotated_cobinding_status/suppressed_merged_S2_to_all_open_and_closed_mnase_peaks_cluster_199_lf_15_rf_15_sec_peak_at_52_lex_300_rex_300_site_peak_1027_4_sam_flag_99~147_cobinding.tsv
    # chr3L:3385499-3385499^199`SRR3133326.3861334_3861334/1_overlapping`99~147#0	MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF............................................................................................................................................FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
    label = line[0:line.find("\t")].split("#")[-1]
    observed_label_list.append(label)
    observed_label_count_dict[label] +=1
    
# randomly shuffle inplace - the observed labels as they are ordered in the file `0`, `1`, 

random.shuffle(observed_label_list)

observed_tf_left =  observed_label_count_dict['3'] + observed_label_count_dict['5']  + observed_label_count_dict['4']  
observed_tf_right =  observed_label_count_dict['1'] + observed_label_count_dict['7']  + observed_label_count_dict['4']  

obsereved_imbalance =  round(abs(math.log (observed_tf_left/(observed_tf_right+1))), 3)
#print("\t".join(["Observed", str(obsereved_imbalance)]))

total_molcules = len(observed_label_list) 

# regenerate original population using bootstraping : use random.choice (with replacement) - with equal proababilities 
n_boots = int(sys.argv[3])

boot_imbalance_list = []
for boot in range(n_boots):
    in_silico_label_count_dict = defaultdict(lambda : 0)
    in_silico_dna_molecules_with_labels = np.random.choice(observed_label_list, size = total_molcules) 
    for k in in_silico_dna_molecules_with_labels: 
        in_silico_label_count_dict[k] +=1
    in_silico_tf_left =  in_silico_label_count_dict['3'] + in_silico_label_count_dict['5']  + in_silico_label_count_dict['4']  
    in_silico_tf_right =  in_silico_label_count_dict['1'] + in_silico_label_count_dict['7']  + in_silico_label_count_dict['4']  
    in_silico_imbalance = abs(math.log((in_silico_tf_left + 1)/(in_silico_tf_right + 1)))
    boot_imbalance_list.append(in_silico_imbalance)
    #print ("\t".join([str(boot), str(in_silico_imbalance)]))


avg_boot = round(np.mean(boot_imbalance_list),3)
std_boot = round(np.std(boot_imbalance_list),3)

t_val = stats.ttest_1samp(boot_imbalance_list, 0)  # 1 is null for balanced
#print("\t".join(["Boot", str(avg_boot), str(std_boot), str(t_val.pvalue)]))

to_write = "\t".join([sys.argv[1], sys.argv[2], str(obsereved_imbalance), str(avg_boot), str(std_boot), str(t_val.pvalue)]) 
out_fp.write(to_write + "\n") 
all_boot_values = "\n".join(map(str, boot_imbalance_list))
out_fp_verb.write(all_boot_values + "\n")


inp_fp_83.close()
inp_fp_99.close()
out_fp.close()
out_fp_verb.close()
