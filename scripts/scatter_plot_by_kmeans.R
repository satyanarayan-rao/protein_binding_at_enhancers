library(ggplot2)
library(ggthemes)
library(grid)
library(Cairo)
args = commandArgs(trailingOnly =  T)
dt = read.table(args[1], header = F, comment.char = "%", stringsAsFactors = F)

dt$lab = unlist(lapply(dt$V1, function(x){
  idx = tail (unlist(strsplit(x, split = "#")),1)
  return (paste0("Cl", " ",  idx))
}  ))
Cairo::CairoPNG(args[2], height = 4, width = 5, units = "in", res = 150)
plt = ggplot(dt, aes(x = V2, y = V3, color = lab)) + 
      geom_point(alpha = 0.3, fill = "white" ) + 
      geom_jitter() + theme_few() + 
      ylab ("Fragment Length [bp]") + xlab("%Orange") + 
  theme(legend.position = "bottom", legend.title = element_blank())
print(plt)
dev.off()

############### Do the histogram of fragment lengths ############ 

dt$reads = unlist(lapply(dt$V1, function(x){
  idx = tail (unlist(strsplit(x, split = "#"))[1])
  return (idx)
}  ))

total_unique_reads = length(unique(dt$reads))
grob_read_count = grobTree(textGrob(paste0("Total Reads: ", total_unique_reads),
                                    x = 0.85, y = 0.95, hjust = 0.5))
Cairo::CairoPNG(args[3], height = 4, width = 5, units = "in", res = 150)

plt = ggplot(dt, aes(x = V3, fill = lab)) + 
      geom_histogram(position = "identity", bins = 25, alpha = 0.5) + 
      xlab("Footprint size [bp]") + ylab ("Count") + 
      geom_rangeframe() + theme_few() + 
      theme(legend.position = "bottom", legend.title = element_blank()) + 
      annotation_custom(grob_read_count)
print (plt)
dev.off()

################# Do the histogram of %Orange ############

Cairo::CairoPNG(args[4], height = 4, width = 5, units = "in", res = 150)
plt = ggplot(dt, aes(x = V2, fill = lab)) + 
      geom_histogram(position = "identity", bins = 25, alpha = 0.5) + 
      xlab("%Orange") + ylab ("Count") + 
      geom_rangeframe() + theme_few() + 
      theme(legend.position = "bottom", legend.title = element_blank()) + 
      annotation_custom(grob_read_count)
print (plt)
dev.off()



