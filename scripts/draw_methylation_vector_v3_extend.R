library(plotrix)
library(stringr)
library(ggplot2)
library(ggthemes)
library(zoo)
source("scripts/methylation_string_to_numeric_vec.R")
source("scripts/matrix_layout.R")
args = commandArgs(trailingOnly = T)
# args[1]: input file (`label`\t`methylation_vec`)
# args[2]: output pdf file
# args[3]: left flank (extended)
# args[4]: right flank (extended)
# args[5]: plot title
# args[6]: plot ylabel
# args[7]: average methylation plot
# args[8]: file to save plot tsv data
# args[9]: strand
# args[10]: MNase data to plot
# args[11]: vertical dash line left (left edge)
# args[12]: vertical dash line right (right edge)
# args[13]: site id 
# args[14]: sam flag
input_dt = read.table(args[1], header = F, stringsAsFactors = F, sep = "\t", comment.char = "%")
names(input_dt) = c("annotation", "methylation_vec")
mnase_dt = read.table(args[10], header = F, sep = "\t", stringsAsFactors = F)
names(mnase_dt) = c("x_tics", "short", "long")
#print (input_dt)
cluster_label = unlist(lapply(input_dt$annotation, 
                  function (x) {unlist(strsplit(x, split = "#"))[2]} ))

input_dt$cluster_label = as.numeric(cluster_label)

input_dt = input_dt[order(input_dt$cluster_label), ]
#print (head(input_dt))

uniq_cluster_label = unique(input_dt$cluster_label)
size_of_each_cluster = lapply(uniq_cluster_label, 
                         function(x) {
                          sub_df = input_dt[which(input_dt$cluster_label == x), ]
                          return(dim(sub_df)[1])
                         })
names(size_of_each_cluster) = uniq_cluster_label
print(size_of_each_cluster)
max_read_of_all_clusters = max (unlist(size_of_each_cluster))
total_reads = sum(unlist(size_of_each_cluster))
print (total_reads)

frac_vec = unlist(lapply(uniq_cluster_label, function (x) {return (size_of_each_cluster[[as.character(x)]]/total_reads)}))
frac_vec_100 = unlist(lapply(uniq_cluster_label, function (x) {return ((size_of_each_cluster[[as.character(x)]]/total_reads)*100)}))
rounded_frac_vec = c()
print (c("--------------\n", frac_vec_100, "\n-----------"))
for (v in frac_vec_100){
    if (v < 1){
        rounded_frac_vec = c (rounded_frac_vec, 1 ) 
    }else{
        rounded_frac_vec = c (rounded_frac_vec, round(v))
    }
}
print (rounded_frac_vec)
print(uniq_cluster_label)
print(str(uniq_cluster_label))
strand = args[9]
lextend = NULL 
rextend = NULL 
lbar = NULL 
rbar = NULL 

if (identical(strand, "+")){
    lextend = as.numeric(args[3])
    rextend = as.numeric(args[4])
    lbar = as.numeric(args[11])
    rbar = as.numeric(args[12])
}else{
    lextend = as.numeric(args[4])
    rextend = as.numeric(args[3])
    lbar = as.numeric(args[11])
    rbar = as.numeric(args[12])
}
x_tic_start = 0 - as.numeric(lextend)
x_tic_stop = as.numeric(rextend)
total_x = x_tic_stop - x_tic_start 
circle_point_size  = (0.55/151)*total_x # 0.5 for 150, so scale it according to that
mid_x = (x_tic_start + x_tic_stop)/2
x_tic_breaks = seq (x_tic_start, x_tic_stop, 25)
cluster_id = 1 
len_label = length(uniq_cluster_label)
cnt = 1 
print (frac_vec)
lout = matrix_layout(0.15, frac_vec)
pdf (args[2], height = 8, width = 6)
#par(mfrow = c (4,1))
#layout(matrix(c(1,1,1,2,2,2,3,3,3,4,4,4), 12, 1, byrow = T))
layout (lout)
par(mar=c(0.5, 0.5, 0.2, 0.2),
     oma = c(2.5, 2.5, 2.5, 2.5))
# plot MNase here
plot (mnase_dt$x_tics, mnase_dt$short, type = "l", main = "", 
      xlim = c (x_tic_start, x_tic_stop), 
      col = "red", frame.plot = F, xaxt = "n", xlab = "", ylab = "", yaxt = "n")
abline(v = c(0 - lbar - 0.5, rbar + 0.5), 
       col = c ("blue", "blue"),
       lty = c(1,1), lwd = c (1,1))
axis (side = 2, at = round(seq(min(mnase_dt$short), max(mnase_dt$short), length.out = 5),2)) 
for (label in uniq_cluster_label){
    non_orange_count = 0
    orange_count = 0 
    plot(1, type="n", xlab="", ylab = label, frame.plot = F,  
         axes = F, yaxt='n',
         main = "", 
         xaxt = "n", 
         xlim=c(x_tic_start, x_tic_stop) , cex = 2,  
         ylim = c(0, size_of_each_cluster[[as.character(label)]] + 1 ),
         #ylim = c(0, 100),
         new = T)
    abline(v = c(0 - lbar - 0.5, rbar + 0.5), 
       col = c ("blue", "blue"),
       lty = c(1,1), lwd = c (1,1))
    if(cnt == len_label){
        axis(side = 1, las = 2, at = x_tic_breaks)
        mtext(paste0(args[13], " (sam flag ", args[14], ")"), side = 1, cex = 0.75, line = 1.5, at = mid_x)
    } else{
        #axis (side = 1, labels = F, color = "white", tick = F, lwd = 0) 
    }
    sub_df = input_dt[which(input_dt$cluster_label == label), ]
    total_lines = dim(sub_df)[1]
    y_val = 0
    for (m_vec in sub_df$methylation_vec){
        y_val = y_val + 1  # used for co-ordinate of y_val
        kp = unlist(strsplit(m_vec, split = ""))
        #print (kp)
        x_tic = x_tic_start
        for(m_s in kp){
            if (identical(m_s, "M")){
                draw.circle(x_tic, y_val + 0.25 , circle_point_size, border = "#bdbdbd", col = "#d7191c", lwd = 0.2)
            }else if (identical (m_s, ".")){
                #x_tic = x_tic + 1
                #next
                non_orange_count = non_orange_count + 1
                draw.circle(x_tic, y_val + 0.25 , circle_point_size, border = "#bdbdbd", col = "white", lwd = 0.2)
            }else if (identical(m_s, tolower(m_s))){
                draw.circle(x_tic, y_val + 0.25 , circle_point_size, lwd = 0.2, border = "#bdbdbd", col = "#bababa")
                non_orange_count = non_orange_count + 1
            }else{
                draw.circle(x_tic, y_val + 0.25 , circle_point_size, col="black", border = "black", lwd = 0.2)
                orange_count = orange_count + 1 
            }
            x_tic = x_tic + 1
        }
    }
    #text(mid_x, max_read_of_all_clusters + 5, paste0(args[5]," ", "Cluster"," ", cluster_id), cex = 1.5)
    cluster_id = cluster_id + 1
    mtext(paste0(size_of_each_cluster[[as.character(label)]], 
         "\n",  round ( orange_count*100/(orange_count + non_orange_count)) , "%" ), 2, cex = 0.5)
    #text(x_tic_start - 5, total_lines/2, args[6], cex = 1.5, srt = 90)
    #mtext (args[6], side = 2, cex = 1.5, line = -1, adj = 0.1)
    cnt = cnt + 1  
}

dev.off()


### print the average methylation status ### 

cl_id = 1
to_plot = NULL
cluster_str = list()
cluster_str_comma = list()
cluster_vec = c()
cluster_comma_vec = c ()
#print(nrow(input_dt))
for (cl in uniq_cluster_label){
    ttp = t(apply(input_dt[which(input_dt$cluster_label == cl), "methylation_vec", drop = F],
               1, convert_string_to_numeric_vec_zero_for_dot))
    #print (names(ttp))
    #ttp = t(ttp)
    #print (head (ttp))
    #print (ttp) 
    #stop("here")
    #avg_m = (colSums(ttp, na.rm = T)*100)/(nrow(input_dt))
    avg_m = (colMeans(ttp, na.rm = T)*100)
    #avg_m = rollmean(avg_m, k = 10, fill = NA) 
    c_s = paste0("Cl", " ", cl_id, "\n", "n = ",
            dim(input_dt[which(input_dt$cluster_label == cl), "methylation_vec", drop = F])[1])
    c_s_comma = paste0("Cl", " ", cl_id, ",", "n = ",
            dim(input_dt[which(input_dt$cluster_label == cl), "methylation_vec", drop = F])[1])
    cluster_vec = c (cluster_vec, c_s)
    cluster_comma_vec = c (cluster_comma_vec, c_s_comma)
    cluster_str[[cl_id]] = c_s 
    cluster_str_comma[[cl_id]] = c_s_comma
    label_vec = rep (c_s, length(avg_m))
    label_comma_vec = rep (c_s_comma, length(avg_m))
    tmp_df = data.frame(avg_m_val = avg_m, clid = label_vec, clid_comma = label_comma_vec,
                  x_tics = as.character(seq(x_tic_start, x_tic_stop)), 
                  stringsAsFactors = F)
    if (is.null(to_plot)){
      to_plot = tmp_df
    }else{
      to_plot = rbind(to_plot, tmp_df)
    }
    cl_id = cl_id + 1 
}
#print (to_plot)
# prepare xtic breaks, keep the minimum, maxium, zero, 
# half of difference zero min of (max and min) 
x_tic_breaks = as.character(c(x_tic_start,
                as.integer(x_tic_start/2), 
                0, as.integer(0 - x_tic_start/2), 
                x_tic_stop)) 
x_tic_breaks = seq (x_tic_start, x_tic_stop, 25)
x_tic_labels = x_tic_breaks
print (x_tic_labels)
print (x_tic_breaks)
#print (x_tic_breaks)
#print (head(to_plot))
to_plot$clid = factor(to_plot$clid, levels = cluster_vec, ordered = T)
to_plot$grey_box_min_x = as.numeric(to_plot$x_tics) + lextend + 0.5 #as.numeric(args[3]) + 0.5
to_plot$grey_box_max_x = as.numeric(to_plot$x_tics) + lextend + 1.5 #as.numeric(args[3]) + 1.5
to_plot$grey_box_min_y = 0 
to_plot$grey_box_max_y = 100
# default alpha = 0.2 
to_plot$alpha_val = 0
x_tics = as.character(seq(x_tic_start, x_tic_stop - 1))

#for (i in x_tics){
#    avg_m_val_at_x = table(to_plot[which(to_plot$x_tics == i), "avg_m_val"], 
#                           useNA = "always")
#    which_is_max = names(which.max(avg_m_val_at_x))
#    if (is.nan(which_is_max)){
#        
#        to_plot[which(to_plot$x_tics == i), "alpha_val"] = 0
#    }
# 
#} 
to_plot[which(to_plot$avg_m_val == 0), "alpha_val"] = 0 
#print (cluster_str)

#print (to_plot)

tt = to_plot[which(to_plot$clid == cluster_str[[1]]), ]$grey_box_min_x

#print (tt)
#print (head(to_plot))
just_cl1_df = to_plot[which (to_plot$clid == cluster_str[[1]]), ]
pdf(args[7], height = 8, width = 8)

#plt = ggplot(to_plot, aes (x = x_tics, y = avg_m_val, color = clid)) + 
#      geom_line(size = 2)  + geom_rangeframe() + theme_few()
#print(plt)

#plt_one_by_one = ggplot(to_plot, aes(x = x_tics, y = avg_m_val, color = clid)) +  
#      geom_line(size = 2)  + geom_rangeframe() + theme_few()
#plt_one_by_one = plt_one_by_one + facet_grid(rows = vars(clid))

plt_one_by_one = ggplot(to_plot, aes(x = x_tics, y = avg_m_val)) + 
                      geom_bar(stat = "identity") + 
                      ylim (c(0,100)) + ylab (paste0("%Methylation", " ", args[6])) +  
                      facet_grid (rows = vars(clid))
#plt_one_by_one = plt_one_by_one + facet_grid(rows = var(to_plot$clid))               

plt_one_by_one = plt_one_by_one + scale_x_discrete( name = "Distance from the peak [bp]", 
                                    limits = as.character(seq(x_tic_start, x_tic_stop - 1)),
                                    breaks = x_tic_breaks, labels = x_tic_labels) + 
                    geom_rangeframe() + theme_few()  + 
                    theme(axis.text.x = element_text(angle = 90))

#plt_one_by_one = plt_one_by_one + annotate ("rect", xmin = 0.5, xmax = 1.5, ymin = 0, ymax = 100, alpha = 0.1) + 
#                   annotate ("rect", xmin = 1.5, xmax = 2.5, ymin = 0, ymax = 100, alpha = 0.1)
for (i in x_tics){ 
    x_min = just_cl1_df[which (just_cl1_df$x_tics == i), "grey_box_min_x"]
    x_max = just_cl1_df[which (just_cl1_df$x_tics == i), "grey_box_max_x"]
    y_min = just_cl1_df[which (just_cl1_df$x_tics == i), "grey_box_min_y"]
    y_max = just_cl1_df[which (just_cl1_df$x_tics == i), "grey_box_max_y"]
    alpha = just_cl1_df[which (just_cl1_df$x_tics == i), "alpha_val"]
    plt_one_by_one = plt_one_by_one + annotate("rect", 
                                       xmin = x_min, xmax = x_max,
                                       ymin = y_min, ymax = y_max,
                                       alpha = alpha) 
}

plt_one_by_one = plt_one_by_one + ggtitle(paste0(args[5], " ", "(total = ", nrow(input_dt) , ")") ) + 
                      theme(plot.title = element_text(hjust = 0.5)) 
print (plt_one_by_one)
dev.off()

to_plot$clid = NULL 
write.table(to_plot, args[8], col.names = T, row.names = F, sep = "\t", quote = F)
