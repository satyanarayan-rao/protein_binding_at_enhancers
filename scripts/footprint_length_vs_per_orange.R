library(Cairo)
args = commandArgs(trailingOnly = T)
dt = read.table(args[1], sep = "\t", comment.char = "%", header = F)
names(dt) = c("peak_id", "per_orange", "footprint_length")

pdf(args[2], height = 4, width = 5)
plot(x = dt$per_orange, y = dt$footprint_length, main = "",
     xlab = "%Orange", ylab = "Footprint size [bp]")
dev.off()
Cairo::CairoPNG(args[3], height = 4, width = 5, units = "in", res = 150)
plot(x = dt$per_orange, y = dt$footprint_length, main = "",
     xlab = "%Orange", ylab = "Footprint size [bp]")
dev.off()

