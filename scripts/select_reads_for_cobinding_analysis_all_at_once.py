import os
import sys
import pickle
import re
mapped_read_dict = pickle.load(open(sys.argv[1], "rb"))
job_id_pkl = pickle.load(open(sys.argv[2], "rb")) 


inp_fp = None 
#mapped_read_dict = None 
peak_id = None 
s_peak_loc = None 
strand = None 
sam_flag = None 
out_fp = None 
lf = None 
rf = None 
for job_k in job_id_pkl: 
    
    inp_fp = open(job_id_pkl[job_k][0])
    #mapped_read_dict = pickle.load(open(sys.argv[2], "rb"))
    peak_id = job_id_pkl[job_k][2]
    s_peak_loc = int(job_id_pkl[job_k][3])
    strand = job_id_pkl[job_k][4]
    sam_flag = job_id_pkl[job_k][5]
    out_fp = open(job_id_pkl[job_k][6], "w")
    lf = int(job_id_pkl[job_k][7])
    rf = int(job_id_pkl[job_k][8])

    #inp_fp = open(sys.argv[1])
    #mapped_read_dict = pickle.load(open(sys.argv[2], "rb"))
    #peak_id = sys.argv[3]
    #s_peak_loc = int(sys.argv[4])
    #strand = sys.argv[5]
    #sam_flag = sys.argv[6]
    #out_fp = open(sys.argv[7], "w")
    #lf = int(sys.argv[8])
    #rf = int(sys.argv[9])
    for line in inp_fp:
        # example line
        # chr3L:632728-632728^10`SRR3133326.1889480_1889480/1_overlapping`99~147	MMMMMMMMMMMMMMMMMMMM.........................................................................................................................................................................................................................................................................................
        read_id = line[line.find("`") + 1: line.find("\t")] 
        read_start = mapped_read_dict[read_id]["start"]
        read_end = mapped_read_dict[read_id]["end"]
        fp_center =  int(line[line.find(":") + 1: line.find("-")])
        
        if (strand == "+") or (strand == "."):
    
            if (fp_center - lf >= read_start) and (fp_center + s_peak_loc + rf <= read_end):
                out_fp.write(line)
            else: 
                continue 
        elif strand == "-": 
            if (fp_center + lf <= read_end) and (fp_center - s_peak_loc - rf >= read_start):
                out_fp.write(line)
            else: 
                continue 
            
    out_fp.close() 
    print ("Finished job {}".format(job_k)) 
    
