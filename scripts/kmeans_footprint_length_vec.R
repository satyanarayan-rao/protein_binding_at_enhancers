args = commandArgs(trailingOnly = T) 

dt = read.table(args[1], sep = "\t", header = F, stringsAsFactors = F, comment.char = "%") 
row.names (dt) = dt$V1
dt$V1 = NULL
lflank = as.integer(args[4])
rflank = as.integer(args[5])
span_from_center = as.integer(args[6])
nclust = as.integer(args[7])
min_of_flanks = min(lflank, rflank)
peak_center = lflank + 1 
dt_sub = NULL 
print ("#############")
print (c(lflank, rflank, nclust, min_of_flanks, peak_center, span_from_center))
if (peak_center - span_from_center >0){ 
    #dt_sub = dt [, seq(peak_center - span_from_center, peak_center + span_from_center)]
    dt_sub = dt
}else {
    #dt_sub = dt [, seq(peak_center - min_of_flanks, peak_center + min_of_flanks)]
    dt_sub = dt 
}
print(head(dt_sub))
set.seed(37)
cl = kmeans(dt_sub, centers = nclust, nstart = 20, iter.max = 200) 

map_df = data.frame(read_id = row.names(dt), cl_id = cl$cluster - 1)
row.names(dt) = paste(row.names(dt), cl$cluster - 1, sep = "#")


write.table(dt, args[2],  row.names = T, col.names = F, quote = F, sep = "\t")
write.table(map_df, args[3], row.names = F, col.names = F, quote = F, sep = "\t")
save (cl, file = args[8])
