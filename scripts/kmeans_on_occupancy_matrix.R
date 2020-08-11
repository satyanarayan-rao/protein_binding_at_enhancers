args = commandArgs(trailingOnly = T)

dt = read.table(args[1], sep = "\t", header = T, stringsAsFactors = F,
                comment.char = "%")
nclust = as.integer(args[2])
#print (dim(dt))
# header example 
# file_name_without_sam_flag	total_reads	total_naked	total_tf	total_nuc	occup_naked	occup_tf	occup_nuc
dt_sub = dt[, c("occup_naked", "occup_tf", "occup_nuc")] 
#print (dim(dt_sub))
#print (head(dt_sub,100))
set.seed(37)

cl = kmeans(dt_sub, centers = nclust, nstart = 20, iter.max = 200)

# reassing cluster basedo on occup_naked percentage ############
cl_id_df = data.frame(actual = cl$cluster) 
center_df = data.frame(cl$centers)

center_sort_df = data.frame(actual = order(center_df$occup_naked), assigned = seq(nclust)) 
cl_id_df$assigned = unlist(lapply(cl_id_df$actual, function(x){
                      return (center_sort_df[which(center_sort_df$actual == x), "assigned"]) 
                         }))
reassinged_center_df_to_write = data.frame(cl$centers)
reassinged_center_df_to_write = reassinged_center_df_to_write[order(reassinged_center_df_to_write$occup_naked),]
row.names(reassinged_center_df_to_write) = seq(nrow(reassinged_center_df_to_write))

cl$reassigned = reassinged_center_df_to_write
df_with_clust = dt
df_with_clust$cl_id = cl_id_df$assigned

write.table(df_with_clust, file = args[3], row.names = F, col.names = T, 
            quote = F, sep = "\t") 
save(cl, file = args[4]) 
write.table (reassinged_center_df_to_write, file = args[5],
             col.names = T, row.names = F, quote = F, sep = "\t")
write.table(data.frame(table (df_with_clust$cl_id)), file = args[6], row.names = F, col.names = F, sep = "\t", quote = F)


## Save the original:

dt_original = dt 
dt_original$cl_id = cl$cluster 
write.table(dt_original, file = args[7], row.names = F, col.names = T, sep = "\t", quote = F )

# save the original centers
write.table(data.frame(cl$centers), file = args[8], row.names = F , col.names = T, sep = "\t", quote = F)
write.table(data.frame(table(cl$cluster)), file = args[9], row.names = F , col.names = F, sep = "\t", quote = F)


################################################################
#df_with_clust$cl_id = cl$cluster


#write.table(df_with_clust, file = args[3], row.names = F, col.names = T, 
#            quote = F, sep = "\t") 
#save(cl, file = args[4])
#write.table (data.frame(cl$centers), file = args[5],
#             col.names = T, row.names = F, quote = F, sep = "\t")
#write.table(data.frame(table (cl$cluster)), file = args[6], row.names = F, col.names = F, sep = "\t", quote = F)
