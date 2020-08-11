library(ggplot2)
library(ggthemes)
library(Cairo)
args = commandArgs(trailingOnly = T)
dt = read.table(args[1], sep = "\t", header = F, comment.char = "%")
print (dim(dt))
names(dt) = c("peak_id", "per_orange")
#pdf(args[2], height = 4, width = 5)
#plt = ggplot (dt, aes(x = per_orange)) + 
#  geom_histogram(position  = "identity", bins = 25, alpha = 0.5) + 
#  theme_few() + geom_rangeframe()  + xlab ("%Orange")
#  #theme(
#  #      #axis.line = element_line(colour = "black"),
#  #      panel.border = element_blank(),
#  #      panel.grid.major = element_blank(),
#  #      panel.grid.minor = element_blank())
#print (plt)
#dev.off()
# 
#Cairo::CairoPNG(args[3], height = 4, width = 5, units = "in", res = 150)
#plt = ggplot (dt, aes(x = per_orange)) + 
#  geom_histogram(position  = "identity", bins = 25, alpha = 0.5) + 
#  theme_few() + geom_rangeframe() + xlab("%Orange") 
#  #theme(
#  #      #axis.line = element_line(colour = "black"),
#  #      panel.border = element_blank(),
#  #      panel.grid.major = element_blank(),
#  #      panel.grid.minor = element_blank())
#print (plt)
#dev.off()
fixed_breaks = seq(0, 100, 10 )
pdf(args[2], height = 4, width = 5)

hist(dt$per_orange, breaks = fixed_breaks, main = "", 
     ylab = "Count", xlab = "%Orange")
dev.off() 
Cairo::CairoPNG(args[3], height = 4, width = 5, units = "in", res = 150)

hist(dt$per_orange, breaks = fixed_breaks, main = "", 
     ylab = "Count", xlab = "%Orange")
dev.off() 

