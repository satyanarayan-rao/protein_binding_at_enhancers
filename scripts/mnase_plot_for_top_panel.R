args = commandArgs(trailingOnly = T) 
# args[1]: mnase data 
# args[2]: output pdf
print (args)
dt = read.table(args[1], sep = "", header = F, stringsAsFactors = F)
pdf(args[2], height = 4, width = 8)
par(mar = c(5,5,2,5))
with (dt, plot(V1, V2, type ="l", col = "blue", ylab = "E.O.M", xlab = NA, ylim = c(min (dt$V2), max(dt$V2) )))
par(new = T, col.lab = "red")
with (dt, plot(x = V1, y= V3, col = "red", type = "l", axes=F, xlab = NA, ylab = NA))
axis (side = 4, col = "red" )
mtext(side = 4, line = 3, "E.O.M")
dev.off()
