library(ggplot2)
library(ggthemes)
library(Cairo)
library(cowplot)


args = commandArgs(trailingOnly = T)
# args[1]: input tsv short
# args[2]: input tsv long
# args[3]: plot title
# args[4]: output pdf
# args[5]: output png
# args[6]: kclust id
dt_short = read.table(args[1], sep = "\t", header = F, stringsAsFactors = F)
dt_long = read.table(args[2], sep = "\t", header = F, stringsAsFactors = F)
cluster_id = paste0("Cluster", args[6])
dt_sub_short =  dt_short[which (dt_short$V3 == cluster_id), ] 
print (head(dt_sub_short))
dt_sub_short["len_type"] = "Short"

dt_sub_long =  dt_long[which (dt_long$V3 == cluster_id), ] 
print (head(dt_sub_long))
dt_sub_long["len_type"] = "Long" 

merged_df = rbind(dt_sub_short, dt_sub_long)

plt = ggplot(merged_df, aes(x = V1, y = V2, color = len_type)) + geom_line() + 
      theme_few() + geom_rangeframe() + xlab ("Distance from peak [bp]") +
      ylab("E.O.M")  + theme(legend.title = element_blank()) + 
      ggtitle(args[3]) + theme(plot.title = element_text(hjust = 0.5))

pdf (args[4], width = 8, height = 4)
print (plt)
dev.off()

Cairo::CairoPNG(args[5], width = 8, height = 4, units = "in", res = 150)
print (plt)
dev.off()
