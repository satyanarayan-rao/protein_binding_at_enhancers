library(plotrix)
library(stringr)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(zoo)
library(Cairo)
library(cowplot)
source("scripts/methylation_string_to_numeric_vec.R")
args = commandArgs(trailingOnly = T)
# args[1]: avg methylation tsv
# args[2]: mnase data
# args[3]: output plot pdf
# args[4]: output plot png
# args[5]: left x-cordinates for grey boxes for example:  "-26,24"; note that `-` has meaning as negative
# args[6]: right x-cordinates for grey boxes for example:  "13,62";  
# for above example two boxes will be created with [xmin1,xmax1] = [-26,13] and [xmin2,xmax2] = [24,62]

grey_box_x_mins = as.numeric(unlist (strsplit(args[5], split = ",")))
grey_box_x_maxs = as.numeric(unlist (strsplit(args[6], split = ",")))
total_boxes = length(grey_box_x_mins) 
dt_methylation = read.table(args[1], header = T, sep = "\t", stringsAsFactors = F)
dt_mnase = read.table(args[2], header = F, sep = "\t", stringsAsFactors = F)
dt_methylation_sub = dt_methylation[, c("x_tics", "avg_m_val", "clid_comma")] 
dt_methylation_sub$clid = unlist(lapply(dt_methylation_sub$clid_comma, function (x){
                              tt = str_replace(x, ",", " (")
                              tt1 = paste0(tt,")")
                              return(tt1) } ))

dt_methylation_sub$clid_comma = NULL 
short_mnase = dt_mnase[, c("V1", "V2")]
names(short_mnase) = c("x_tics", "avg_m_val")
short_mnase$clid = "MNase (Short)"

to_plot_df = rbind(short_mnase, dt_methylation_sub)
individual_clusters = unique(to_plot_df$clid) 
plt_list = list() 
cnt = 1
total_sub_plots = length(individual_clusters)
for (cl_id in individual_clusters){
     plt = NULL
     tmp_df = to_plot_df[which(to_plot_df$clid ==cl_id), ]
     if (identical (cl_id, "MNase (Short)")){
         if (cnt !=total_sub_plots){
             plt = ggplot(tmp_df, aes(x = x_tics, y = avg_m_val)) + 
                   geom_line( col = "red") + theme_few() + 
                   theme(axis.title.x = element_blank(), 
                         axis.ticks.x=element_blank()) + ylab(cl_id) + 
                   theme(axis.text.x = element_blank()) + 
                   geom_area(fill = "red", alpha = 0.1) + 
                   theme(plot.margin = unit(c(0.15, 0, 0, 0), "cm")) 
             for (gbox in seq (total_boxes)){
                   plt = plt + annotate ("rect", xmin = grey_box_x_mins[gbox], 
                             xmax = grey_box_x_maxs[gbox], 
                             ymin = -Inf, ymax = Inf, alpha = 0.2, border = NA)  
             }
         }else{
             plt = ggplot(tmp_df, aes(x = x_tics, y = avg_m_val)) + 
                   geom_line( col = "red") + theme_few() + xlab("") + ylab(cl_id) +
                   theme(axis.text.x = element_blank()) + 
                   theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) 
             for (gbox in seq (total_boxes)){
                   plt = plt + annotate ("rect", xmin = grey_box_x_mins[gbox], 
                             xmax = grey_box_x_maxs[gbox], 
                             ymin = -Inf, ymax = Inf, alpha = 0.2, border = NA)  
             }

         }
         
     }else{
         if (cnt!=total_sub_plots){
         plt = ggplot(tmp_df, aes(x = x_tics, y = avg_m_val)) + 
               geom_bar(stat = "identity", color = "orange") + geom_rangeframe () +
               theme_few() + xlab("") + ylab(cl_id) + 
               theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
                     axis.ticks.x=element_blank()) + 
               theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + ylim (c(0,100))
         for (gbox in seq (total_boxes)){
               plt = plt + annotate ("rect", xmin = grey_box_x_mins[gbox], 
                         xmax = grey_box_x_maxs[gbox], 
                         ymin = -Inf, ymax = Inf, alpha = 0.2, border = NA)  
             }
              
         }else{
             plt = ggplot(tmp_df, aes(x = x_tics, y = avg_m_val)) + 
                   geom_bar(stat = "identity", color = "orange") + 
                   geom_rangeframe () + theme_few() + ylab(cl_id) + 
                   theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + 
                   xlab("Distacne from peak (0) [bp]") + ylim(c(0,100)) 
             for (gbox in seq (total_boxes)){
                   plt = plt + annotate ("rect", xmin = grey_box_x_mins[gbox], 
                             xmax = grey_box_x_maxs[gbox], 
                             ymin = -Inf, ymax = Inf, alpha = 0.2, border = NA)  
             }
                   
         }
     }
     
     plt_list[[cnt]] = plt 
     cnt = cnt + 1
}

#plt_all = gridExtra::grid.arrange (plt_list[[1]], plt_list[[2]], plt_list[[3]], plt_list[[4]], plt_list[[5]], plt_list[[6]], plt_list[[7]])
#print (tail(to_plot_df))
gl = lapply(plt_list, ggplotGrob)
pdf(args[3], width = 3, height = 9)
#plt_one_by_one = ggplot(to_plot_df, aes(x = x_tics, y = avg_m_val)) + 
#                      geom_bar(stat = "identity") + 
#                      ylab (paste0("%Fooprint", " ", args[6])) + 
#                      facet_grid (rows = vars(clid)) + geom_rangeframe() + theme_few()
#print (plt_one_by_one)
#gridExtra::grid.arrange (plt_list[[1]], plt_list[[2]], plt_list[[3]], plt_list[[4]], plt_list[[5]], plt_list[[6]], plt_list[[7]], ncol = 1) 
plot_grid(plotlist = gl, ncol  = 1, align = 'v')
dev.off()

Cairo::CairoPNG(args[4], width = 3, height = 9, units = "in", res = 150)
#gridExtra::grid.arrange (plt_list[[1]], plt_list[[2]], plt_list[[3]], plt_list[[4]], plt_list[[5]], plt_list[[6]], plt_list[[7]], ncol = 1) 
plot_grid(plotlist = gl, ncol  = 1, align = 'v')
dev.off()

