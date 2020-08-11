import sys
inp_fp = open(sys.argv[1])
out_fp = open(sys.argv[2], "w")

for line in inp_fp:
    o_vec = [] 
    d_loc = line.find("\t") 
    read_id = line[0:d_loc]
    vec = line[d_loc + 1: len(line)]
    for c in range (len(vec) - 1): 
        if vec[c] ==  "." : 
            o_vec.append("0") 
        elif vec[c] == "M":
            o_vec.append("-1")
        elif vec[c] == vec[c].lower():
            o_vec.append("1")
        elif vec[c] == vec[c].upper():
            o_vec.append("2")
            
        
    to_write = read_id + "\t" + "\t".join(o_vec)
    out_fp.write(to_write + "\n")

out_fp.close() 
inp_fp.close()     
    
