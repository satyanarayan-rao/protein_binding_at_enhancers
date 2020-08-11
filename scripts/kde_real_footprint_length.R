args = commandArgs(trailingOnly = T) 
options(bitmapType='cairo')
dt = read.table(args[1], sep = "", header = T, stringsAsFactor = F)
dt_colnames = names(dt)
col_to_group = as.numeric(args[2])
col_to_plot = as.numeric(args[3])
bw = as.numeric(args[4])
color_data = read.table(args[5], header = F, stringsAsFactor = F)
label_prefix = args[6]
x_label = args[7]
y_label = args[8]
line_styles = read.table(args[9], header = F, stringsAsFactor = F)
nclust = as.numeric(args[12])
dt_sub = dt [, c(col_to_group, col_to_plot)]
print (head(dt_sub))

dt_sub_l = lapply(unique(dt_sub[[dt_colnames[col_to_group]]]), 
                  function (x) { 
		      dt_sub[which(dt_sub[[dt_colnames[col_to_group] ]]==x), dt_colnames[col_to_plot] ] })
names(dt_sub_l) = unique(dt_sub[[dt_colnames[col_to_group]]])

hist_x = list() 
hist_y = list()
line_labels = list() 
line_labels_vec = c()
color_vals = c()
line_vals = c()
min_y = 100
max_y = -1
min_x = 10000
max_x = -1 
counter = 1
to_save_df = NULL
for(l in names(dt_sub_l)){
    lab = paste0(label_prefix, " ", l)
    line_labels[[l]] =  lab
    line_labels_vec = c (line_labels_vec, lab)
    color_vals = c (color_vals, color_data[counter, "V1"])
	  line_vals = c (line_vals, line_styles[counter, "V1"])
    #print(dt_sub_l[[l]])
    hist_group = density(dt_sub_l[[l]], bw = 20, from = 0, to = 200, n = 30)
    tmp_df = data.frame(x_val = hist_group$x, y_val = hist_group$y, 
                        cl_id = rep(l, length(hist_group$x)))
    if (is.null(to_save_df)){
        to_save_df = tmp_df
    }else{
        to_save_df = rbind(to_save_df, tmp_df)
    }
    #hist_group = density(dt_sub_l[[l]],  n = 30)
    hist_x[[l]] = hist_group$x 
    hist_y[[l]] = hist_group$y
		if (max(hist_y[[l]]) > max_y ){
		    max_y = max(hist_y[[l]])
		}
		if (min (hist_y[[l]]) < min_y) {
		    min_y = min(hist_y[[l]])
		} 
    if (max(hist_x[[l]]) > max_x ){
        max_x = max(hist_x[[l]]) 
    }
    if (min(hist_x[[l]]) < min_x){
       min_x = min(hist_x[[l]]) 
    }
    counter = counter + 1

}
write.table(to_save_df, args[15], row.names = F, col.names = F, quote = F, sep = "\t")
print (hist_y)
counter = 1

matrix_nclust = NULL
if (nclust%%2){
    matrix_nclust = nclust + 1 
} else{
    matrix_nclust = nclust
}


png(args[10], height = 3*(matrix_nclust/2), width = 6, units = "in", res = 300)

#par (mfrow =c(3,2), byrow = T)
layout(matrix(seq(matrix_nclust), matrix_nclust/2, 2, byrow = T ))
for (l in names(dt_sub_l)){
  plot(x = hist_x[[l]], y = hist_y[[l]], 
       type = 'l', col = color_data[counter, "V1"],
       xlab = x_label, ylab = y_label, xlim = c(min_x, max_x), lty = line_styles[counter, "V1"],
       main = paste0(label_prefix, counter))
  color_vals = c (color_vals, color_data[counter, "V1"])
	line_vals = c (line_vals, line_styles[counter, "V1"])
  counter = counter + 1 
}

dev.off()
pdf(args[11], height = 3*(matrix_nclust/2), width = 6)
for (l in names(dt_sub_l)){
  plot(x = hist_x[[l]], y = hist_y[[l]], 
       type = 'l', col = color_data[counter, "V1"],
       xlab = x_label, ylab = y_label, xlim = c(min_x, max_x), lty = line_styles[counter, "V1"],
       main = paste0(label_prefix, counter))
  color_vals = c (color_vals, color_data[counter, "V1"])
	line_vals = c (line_vals, line_styles[counter, "V1"])
  counter = counter + 1 
}

dev.off()

##################### all kde in one plot ############################# 
png (args[13], width = 7.5, height = 5, units = "in", res = 300)
counter = 1
color_vals = c()
line_vals = c() 
for (l in names(dt_sub_l)){
    if (counter == 1){
        plot(x = hist_x[[l]], y = hist_y[[l]], 
             type = 'l', col = color_data[counter, "V1"],
             xlab = x_label, ylab = y_label, ylim = c(min_y, max_y), xlim = c(min_x, max_x), lty = line_styles[counter, "V1"])
        color_vals = c (color_vals, color_data[counter, "V1"])
	line_vals = c (line_vals, line_styles[counter, "V1"])
    }else{
        lines(x = hist_x[[l]], y = hist_y[[l]], col = color_data[counter, "V1"], lty = line_styles[counter, "V1"])
        color_vals = c (color_vals, color_data[counter, "V1"])
	line_vals = c (line_vals, line_styles[counter, "V1"])
    }
    counter = counter + 1 
}

legend ("topright", legend = line_labels_vec, col = color_vals, lty = line_vals)
dev.off()


################ pdf  ########### 
pdf (args[14], width = 7.5, height = 5)
counter = 1
color_vals = c()
line_vals = c() 
for (l in names(dt_sub_l)){
    if (counter == 1){
        plot(x = hist_x[[l]], y = hist_y[[l]], 
             type = 'l', col = color_data[counter, "V1"],
             xlab = x_label, ylab = y_label, ylim = c(min_y, max_y), xlim = c(min_x, max_x), lty = line_styles[counter, "V1"])
        color_vals = c (color_vals, color_data[counter, "V1"])
	line_vals = c (line_vals, line_styles[counter, "V1"])
    }else{
        lines(x = hist_x[[l]], y = hist_y[[l]], col = color_data[counter, "V1"], lty = line_styles[counter, "V1"])
        color_vals = c (color_vals, color_data[counter, "V1"])
	line_vals = c (line_vals, line_styles[counter, "V1"])
    }
    counter = counter + 1 
}

legend ("topright", legend = line_labels_vec, col = color_vals, lty = line_vals)
dev.off()
################# all footprint length in one ##############
fixed_breaks = seq(0, 300, length.out = 100)
h_total = hist(dt[, c(col_to_plot)], breaks = fixed_breaks, plot = F)
pdf (args[16], width = 5, height = 4)
hist(dt[, c(col_to_plot)], breaks = fixed_breaks, main = "", 
               ylab = paste0("Count in ", args[18]), xlab = "Footprint size [bp]")
dev.off()
png(args[17], width = 5, height = 4, units = "in", res = 150)
hist(dt[, c(col_to_plot)], breaks = fixed_breaks, main = args[18], 
               ylab = paste0( "Count in ", args[18]), xlab = "Footprint size [bp] ")
dev.off()
################## cap 200 ####################################
fixed_breaks = seq(0, 200, length.out = 70)
to_plot = dt[, c(col_to_plot)]
gt_200 = which(to_plot > 200)
to_plot[gt_200] = 200
pdf (args[19], width = 5, height = 4)
hist(to_plot, breaks = fixed_breaks, main = "", 
               ylab = paste0("Count in ", args[18]), xlab = "Footprint size [bp]", 
               xaxt = "n") 
where = seq(0,200, 20) 
txt_where = as.character(where)
txt_where[length(txt_where)] = paste0("\u2265 ", txt_where[length(txt_where)])
axis(1, at = where, labels = txt_where, las = 2)
dev.off()

png(args[20], width = 5, height = 4, units = "in", res = 150)
hist(to_plot, breaks = fixed_breaks, main = "", 
               ylab = paste0("Count in ", args[18]), xlab = "Footprint size [bp]", 
               xaxt = "n") 
where = seq(0,200, 20) 
txt_where = as.character(where)
txt_where[length(txt_where)] = paste0("\u2265 ", txt_where[length(txt_where)])
axis(1, at = where, labels = txt_where, las = 2)
dev.off()
