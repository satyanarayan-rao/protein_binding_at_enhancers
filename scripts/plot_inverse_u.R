library(Cairo)
library(ggplot2)
library(ggthemes)
args = commandArgs(trailingOnly = T) 

dt = read.table(args[1], sep = "", header = F, comment.char = "%") 

row.names(dt) = dt$V1 
dt$V1 = NULL
dt_copy = dt 
dt_copy$cl_id = unlist(lapply (row.names(dt), function(x){
           tail (unlist (strsplit(x, split = "#")),1) })) 

all_cl_ids = sort(as.numeric(unique(dt_copy$cl_id))) # making it first integer and then sorting
all_cl_ids = as.character(all_cl_ids)
remap_cl = list() 
cl_vec = c()
cnt = 1
for (i in all_cl_ids){
     remap_cl[[i]] = cnt 
     cl_vec = c(cl_vec, cnt)
     cnt = cnt + 1 
     
}
cl_vec = as.character(cl_vec)
print (cl_vec)
print (remap_cl)
dt_copy$remapped_cl = unlist(lapply(dt_copy$cl_id, function (x){
            return (remap_cl[[x]])
                 })) 
print (head(dt_copy))
# calculate colMeans of df  
mean_list  = list() 
median_list = list() 
cluster_count_list = list()
for (cl in cl_vec){
    dt_sub = dt [ which(dt_copy$remapped_cl == cl),  ] 
    mean_list[[cl]] = colMeans(dt_sub) 
    median_list[[cl]] = apply (dt_sub, 2, median)  
    cluster_count_list[[cl]] = nrow(dt_sub)
}

nclust = as.integer(args[5])
matrix_nclust = NULL
if (nclust%%2){
    matrix_nclust = nclust + 1 
} else{
    matrix_nclust = nclust
}
x_vals = seq (0 - as.numeric(args[3]), as.numeric (args[4])) 

pdf(args[2], height = 3*nclust, width = 6)
layout(matrix(seq(nclust*2), nclust, 2, byrow = T ))
#par(mar=c(0.5, 0.5, 0.2, 0.2),
#     oma = c(2.5, 2.5, 2.5, 2.5))
for (cl in cl_vec){
    plot(x = x_vals, y = mean_list[[cl]],
         type = "l", lwd = 2, col = "red",
         main = paste0("Cluster ", cl, " (n = ", cluster_count_list[[cl]], ")"),
         xlab = "MNase boundary", ylab = "Mean footprint length")
    plot(x = x_vals, y = median_list[[cl]],
         type = "l", lwd = 2, col = "red",
         main = paste0("Cluster ", cl, " (n = ", cluster_count_list[[cl]], ")"),
         xlab = "MNase boundary", ylab = "Median footprint length")
}
dev.off ()
Cairo::CairoPNG(args[6], height = 3*nclust, width = 6, units = "in", res = 150)
layout(matrix(seq(nclust*2), nclust, 2, byrow = T ))
#par(mar=c(0.5, 0.5, 0.2, 0.2),
#     oma = c(2.5, 2.5, 2.5, 2.5))
for (cl in cl_vec){
    plot(x = x_vals, y = mean_list[[cl]],
         type = "l", lwd = 2, col = "red",
         main = paste0("Cluster ", cl, " (n = ", cluster_count_list[[cl]], ")"),
         xlab = "MNase boundary", ylab = "Mean footprint length")
    plot(x = x_vals, y = median_list[[cl]],
         type = "l", lwd = 2, col = "red",
         main = paste0("Cluster ", cl, " (n = ", cluster_count_list[[cl]], ")"),
         xlab = "MNase boundary", ylab = "Median footprint length")
}
dev.off ()
