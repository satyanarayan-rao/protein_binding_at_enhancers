import os
import sys
import Pycluster as PC
import nltk
import random
# sys.argv[1]: input methylation vector file
# sys.argv[2]: site to select for
# sys.argv[3]: output methylation vector file with sam_flag 83~163 file 
# sys.argv[4]: output methylation vector file with sam_flag 99~147 file
# sys.argv[5]: number of cluster
to_write_sam_flag_83_163 = open (sys.argv[3], "w")
to_write_sam_flag_99_147 = open (sys.argv[4], "w")
site_reads_flag_83_163 = []
site_reads_flag_99_147 = []
site_reads_flag_83_163_label = []
site_reads_flag_99_147_label = []
input_fp = open (sys.argv[1])
random.seed(71)
for line in input_fp: 
    if line.find(sys.argv[2]) != -1: 
        if line.find("83~163") != -1: 
            #to_write_sam_flag_83_163.write(line)
            tab_loc = line.find("\t")
            site_reads_flag_83_163.append(line[tab_loc + 1: len(line) - 1])
            site_reads_flag_83_163_label.append(line[0:tab_loc])
        elif line.find("99~147") != -1: 
            tab_loc = line.find("\t")
            site_reads_flag_99_147.append(line[tab_loc + 1: len(line) - 1])
            site_reads_flag_99_147_label.append(line[0:tab_loc])
    else:
        continue
print (len(site_reads_flag_83_163))
dist_83_163 = [nltk.distance.edit_distance(site_reads_flag_83_163[i],
                                  site_reads_flag_83_163[j])
               for i in range(0, len(site_reads_flag_83_163))
               for j in range(0, i)]
dist_99_147 = [nltk.distance.edit_distance(site_reads_flag_99_147[i],
                                  site_reads_flag_99_147[j])
               for i in range(0, len(site_reads_flag_99_147))
               for j in range(0, i)]
print (len(site_reads_flag_99_147))

nclust = int(sys.argv[5])
labels_83_163, error_83_163, nfound_83_163 = PC.kmedoids(dist_83_163, nclusters = nclust)
labels_99_147, error_99_147, nfound_99_147 = PC.kmedoids(dist_99_147, nclusters = nclust)
cluster_83_163 = dict()
for vec, cl_label, vec_label in zip (site_reads_flag_83_163, 
                                     labels_83_163,
                                     site_reads_flag_83_163_label): 
    read_with_label = vec + "@" + vec_label
    cluster_83_163.setdefault(cl_label, []).append(read_with_label) 
    to_write = vec_label + "#" + str(cl_label) + "\t" + vec
    to_write_sam_flag_83_163.write(to_write + "\n")
cluster_99_147 = dict()
for vec, cl_label, vec_label in zip (site_reads_flag_99_147, 
                                     labels_99_147,
                                     site_reads_flag_99_147_label): 
    read_with_label = vec + "@" + vec_label
    cluster_99_147.setdefault(cl_label, []).append(read_with_label) 
    to_write = vec_label + "#" + str(cl_label) + "\t" + vec
    to_write_sam_flag_99_147.write(to_write + "\n")


to_write_sam_flag_83_163.close()
to_write_sam_flag_99_147.close()
