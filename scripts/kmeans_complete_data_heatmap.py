import sys
import subprocess

tsv_file=sys.argv[1]
basic_commands_fp = open (sys.argv[2])
out_fp = open (sys.argv[3], "w")
tsv_fp = open (tsv_file)

cnt_fp = open (sys.argv[6])
cnt_dict = {} 
total_clusters = 0
cl_line_dict = {} 
cumsum_dict = {} 
for line in cnt_fp:
    total_clusters  +=1
    cl_id_and_cnt = line.strip().split()
    #cnt_dict[cl_id_and_cnt[0]] = '"' + str (total_clusters) + ' (n = ' +  "%4s"%(cl_id_and_cnt[1]) + ')' + '"' 
    cnt_dict[cl_id_and_cnt[0]] =  str (total_clusters) + ' (n = ' +  "%4s"%(cl_id_and_cnt[1]) + ')' 
    cl_line_dict[cl_id_and_cnt[0]] = cl_id_and_cnt[1]

for cl in range(2, total_clusters + 2):
    csum = 0 
    for c in range(1, cl):
        csum += int(cl_line_dict[str(c)])
    cumsum_dict[cl - 1] = csum

#line_breaks_file=sys.argv[2]
out = subprocess.Popen(['wc', '-l', tsv_file], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)  
stdout, stderr = out.communicate()
num_lines = int (stdout.decode().split()[0]) - 1  
line_headers = tsv_fp.readline().strip().split()[5:8]
len_line_headers = len(line_headers)

# what we want is the string of this format: ("x" 0, "y" 1)
line_headers_with_quotes = [ '"' + i + '"' for i in line_headers]
gnuplot_xlab_string = [ j + " " + str(i) for i,j in zip (list(range(len_line_headers)), line_headers_with_quotes )]
gnuplot_ylab_string = [ cnt_dict[str(i)] + " " + str (cumsum_dict[i]/2) for i in range (1, total_clusters + 1)]

xtics_str = "(" + ",".join(gnuplot_xlab_string) + ")" + " " + "rotate by -45"
ytics_str = "(" + ",".join(gnuplot_ylab_string) + ")"

       
  
for line in basic_commands_fp: 
    out_fp.write(line)

#out_fp.write("# ---------------- writing values in cell --------------#\n")
#tsv_fp.seek(0)
#h = tsv_fp.readline() # ignore the header
#row_id = 1
#label_id = 1 
#for line in tsv_fp:
#    col_id = 0
#    lv = line.strip().split()
#    for v in lv:
#        rv = round (float(v), 2)
#        out_fp.write("set label {l_id} at {x}, {y} '{text}' center front\n".format (
#           l_id = label_id, x = col_id, y = row_id, text = rv)) 
#        col_id +=1 
#        label_id +=1 
#    row_id +=1 

out_fp.write ("set output '{f}'\n".format(f=sys.argv[4]))
out_fp.write ("set yrange [0.5:{y}]\n".format(y = num_lines + 0.5)) 
#out_fp.write("set ytics 1,1,{nlines}\n".format(nlines = num_lines))
#out_fp.write("set ytics {labels}\n".format(labels = ytics_str))
out_fp.write ("set xtics {labels}\n".format (labels = xtics_str))


# write line numbers

for i in range(1, total_clusters):
    line_loc = 0 
    for j in range(i):
        line_loc += int(cl_line_dict[str(j + 1 )]) 
    
    out_fp.write("set arrow {l_id} from -2.5, {line_loc} to {nclust}, {line_loc} nohead lc rgb '#000000' front\n".format(l_id = i, line_loc = line_loc, nclust=total_clusters - 0.5)) 
out_fp.write("set arrow {l_id} from -2.5, {line_loc} to {nclust}, {line_loc} nohead lc rgb '#000000' back\n".format(l_id = total_clusters, line_loc = num_lines, nclust = total_clusters - 0.5))
for cl in range(1, total_clusters + 1):
    if cl == 1:
        out_fp.write("set label {l_id} '{lab_str}'  at -1.5, {ypos} right front\n".format(l_id = cl, lab_str = cnt_dict[str(cl)],  ypos = cumsum_dict[cl]/2)) 
    else:
        out_fp.write("set label {l_id} '{lab_str}'  at -1.5, {ypos} right front\n".format(l_id = cl, lab_str = cnt_dict[str(cl)],  ypos = (cumsum_dict[cl] + cumsum_dict[cl - 1])/2)) 
   

out_fp.write ("set title '{f}'\n".format(f=sys.argv[5]))
out_fp.write ("plot '<cut -f2-9 {f} | sort -k8,8n | cut -f5-7' matrix with image notitle\n".format(f=sys.argv[1]))

basic_commands_fp.close() 
out_fp.close()
