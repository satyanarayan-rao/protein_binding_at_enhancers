args = commandArgs(trailingOnly = T) 
options(bitmapType='cairo')
# Rscript distplot_based_on_one_column_with_header.R <input_file> <columd id to use for group plot>  <column id to use for distribution plotting> <number of breaks> <color seletion file> <label prefix> <x label text> <y label text>
# Rscript  $NGS_SCRIPTS_DIR/distplot_based_on_one_column_with_header.R motif2cl_id.tsv_with_patch_dist.tsv 2 3 20 $NGS_SCRIPTS_DIR/r_colors.list Cluster "Distance from motif center (bp)" "Frequency"
dt = read.table(args[1], sep = "", header = T, stringsAsFactor = F)
print (args)
dt_colnames = names(dt)
col_to_group = as.numeric(args[2])
col_to_plot = as.numeric(args[3])
num_breaks = as.numeric(args[4])
color_data = read.table(args[5], header = F, stringsAsFactor = F)
label_prefix = args[6]
x_label = args[7]
y_label = args[8]
line_styles = read.table(args[9], header = F, stringsAsFactor = F)

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
min_y = 100
max_y = -1
min_x = 200
max_x = -1 
for(l in names(dt_sub_l)){
    lab = paste0(label_prefix, " ", l)
    line_labels[[l]] =  lab
    line_labels_vec = c (line_labels_vec, lab)
    hist_group = hist(dt_sub_l[[l]], breaks = num_breaks, plot = F)
    hist_x[[l]] = hist_group$mids 
    hist_y[[l]] = hist_group$density
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

}

png(paste0(args[1], "-dist.png"), height = 4, width = 6, units = "in", res = 300)
#pdf(paste0(args[1], "-dist.pdf"), height = 6, width = 9)

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
counter = 1 
png(paste0(args[1], "-dist-ind.png"), height = 9, width = 6, units = "in", res = 300)

#par (mfrow =c(3,2), byrow = T)
layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow = T ))
for (l in names(dt_sub_l)){
  plot(x = hist_x[[l]], y = hist_y[[l]], 
       type = 'l', col = color_data[counter, "V1"],
       xlab = x_label, ylab = y_label, xlim = c(min_x, max_x), lty = line_styles[counter, "V1"],
       main = paste0("Cl ", counter))
  color_vals = c (color_vals, color_data[counter, "V1"])
	line_vals = c (line_vals, line_styles[counter, "V1"])
  counter = counter + 1 
}

dev.off()

############### kde based 
print ("hellooo...\n")
hist_x = list() 
hist_y = list()
min_y = 100
max_y = -1
min_x = 10000
max_x = -1 
for(l in names(dt_sub_l)){
    lab = paste0(label_prefix, " ", l)
    line_labels[[l]] =  lab
    line_labels_vec = c (line_labels_vec, lab)
    print(dt_sub_l[[l]])
    hist_group = density(dt_sub_l[[l]], bw = 20, from = 0, to = 200, n = 30)
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

}
print (hist_y)
counter = 1
png(paste0(args[1], "-kde-ind.png"), height = 9, width = 6, units = "in", res = 300)

#par (mfrow =c(3,2), byrow = T)
layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow = T ))
for (l in names(dt_sub_l)){
  plot(x = hist_x[[l]], y = hist_y[[l]], 
       type = 'l', col = color_data[counter, "V1"],
       xlab = x_label, ylab = y_label, xlim = c(min_x, max_x), lty = line_styles[counter, "V1"],
       main = paste0("Cl ", counter))
  color_vals = c (color_vals, color_data[counter, "V1"])
	line_vals = c (line_vals, line_styles[counter, "V1"])
  counter = counter + 1 
}

dev.off()

