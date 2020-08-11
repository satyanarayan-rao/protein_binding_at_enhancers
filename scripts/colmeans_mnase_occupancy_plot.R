library(ggplot2)
library(ggthemes)
library(Cairo)

args = commandArgs(trailingOnly = T)
# args[1]: input tsv
# args[2]: output pdf
# args[3]: output png
dt = read.table(args[1], sep = "\t", header = F, stringsAsFactors = F)

plt = ggplot(dt, aes(x = V1, y = V2, color = V3)) + geom_line() + 
      theme_few() + geom_rangeframe() + xlab ("Distance from peak [bp]") +
      ylab("E.O.M")  + theme(legend.title = element_blank()) + 
      ggtitle(args[4]) + theme(plot.title = element_text(hjust = 0.5))

pdf (args[2], width = 8, height = 4)
print (plt)
dev.off()

Cairo::CairoPNG(args[3], width = 8, height = 4, units = "in", res = 150)
print (plt)
dev.off()
