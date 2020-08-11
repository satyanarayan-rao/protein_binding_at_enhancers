import sys
inp_fp = open(sys.argv[1])
out_fp = open(sys.argv[2], "w")

for line in inp_fp:
    o_vec = [] 
    for c in range (len(line) - 1): 
        if line[c] ==  "." : 
            o_vec.append("0") 
        elif line[c] == "M":
            o_vec.append("-1")
        elif line[c] == "F":
            o_vec.append("1")
        
    to_write = "\t".join(o_vec)
    out_fp.write(to_write + "\n")

out_fp.close() 
inp_fp.close()     
    
