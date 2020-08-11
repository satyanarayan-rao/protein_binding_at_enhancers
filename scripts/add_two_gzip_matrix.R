library(R.utils)
library(stringr)
args = commandArgs(trailingOnly = T)
dat1 = read.table(args[1], sep = "", header = F, stringsAsFactors = F)
dat2 = read.table(args[2], sep = "", header = F, stringsAsFactors = F)

row.names(dat1) = dat1$V1
row.names(dat2) = dat2$V1
dat1[, 1] = NULL
dat2[, 1] = NULL
# do the intersection of row ids
common_rows = intersect (row.names(dat1), row.names(dat2))
summed_df = dat1[common_rows, ] + dat2[common_rows, ]
summed_df = cbind(common_rows = common_rows, summed_df)
unzipped_fname = str_replace(args[3], ".gz$", "")
write.table(summed_df, unzipped_fname, row.names = F, col.names = F, quote = F, sep = "\t")
gzip(unzipped_fname)
