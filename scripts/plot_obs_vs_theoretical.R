library(ggplot2)
library(ggthemes)
library(Cairo)
library(zoo)
args = commandArgs(trailingOnly = T)

# ------------------- #

list_of_files = unlist(strsplit(args[1], split = " "))
list_of_labels = unlist(strsplit(args[2], split = "@"))
out_pdf = args[3]
out_png = args[4] 
to_plot_df = NULL
for (i in seq(length(list_of_files))){
    tmp_df = read.table(list_of_files[i], sep = "\t", header = F, stringsAsFactors = F)
    tmp_df$V4 = list_of_labels[i]  
    tmp_df = tmp_df[tmp_df$V1>0 & tmp_df$V1<=300,] 
    tmp_df = tmp_df[complete.cases(tmp_df), ]
    if (grepl("Obs", list_of_labels[i])){
        tmp_df$rmean = rollmean(tmp_df$V3, k = 5, fill = 0)
        tmp_df = tmp_df[tmp_df$V1>0, ] 
    }
    if (is.null(to_plot_df)){
        to_plot_df = tmp_df
    }else{
        if (grepl("Theoretical", list_of_labels[i])){
        	tmp_df$rmean = rollmean(tmp_df$V3, k = 5, fill = 0)
        }
        to_plot_df = rbind(to_plot_df, tmp_df)
    }
    #print (to_plot_df)
}
to_plot_df$label = paste(to_plot_df$V4, to_plot_df$V5, sep = "-")
to_rename_legned_keys = ("")

plt = ggplot(to_plot_df, aes(x = V1, y = rmean, color = V4)) + geom_line() 

pdf(out_pdf, width = 6, height = 5)
print(plt)
dev.off()

Cairo::CairoPNG(out_png, width = 6, height = 5, units = "in", res = 150)
print(plt)
dev.off()

write.table(to_plot_df, args[5], row.names = F, col.names = F, quote = F, sep = "\t")
