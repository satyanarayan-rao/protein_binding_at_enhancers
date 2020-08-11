library(stringr)
library(ggplot2)
library(ggthemes)
library(Cairo)
# args[1]: input file list ; .Rdata files from kmean 
# args[2]: label list for input files; cluster id associated with .Rdata file
# args[3]: output tsv file ; to store the wss
# args[4]: output png 
options(error=traceback)
args = commandArgs(trailingOnly = T)
file_list = unlist(strsplit(args[1], split = " "))
label_list = unlist(strsplit(args[2], split = "@"))
cnt = 1 
wss_vec = c()
cl_vec = c() 
for (f in file_list){
    load(f)
    total_wss = sum(cl$withinss)
    wss_vec = c(wss_vec, total_wss)
    cl_vec = c(cl_vec, as.integer(label_list[cnt])) 
    cnt = cnt + 1 
}
wss_df = data.frame(cl_id = cl_vec, total_wss = wss_vec) 
write.table(wss_df, file = args[3], 
            sep = "\t", row.names = F, col.names = T, quote = F)
Cairo::CairoPNG(args[4], width = 8, height = 5, units = "in", res = 150)  
plot(y = wss_df$total_wss, x = wss_df$cl_id,  
      main = paste(args[5], "kmeans elbow", sep = " "),
      xlab = "Number of clusters", ylab = "Total WSS", 
      lty = 1, lwd = 1, type = "b", 
      col = "red", pch = 19, xaxt = "n", frame = F)
axis(side = 1, at = wss_df$cl_id, las = 2)
dev.off() 
