library(ggplot2)
library(ggthemes)
library(cowplot)
library(Cairo)
library(grid)

args = commandArgs(trailingOnly = T)

dt = read.table(args[1], sep = "\t", header = F, stringsAsFactors = F)

dt_p = dt[dt$V1=="primary", ]
dt_s = dt[dt$V1=="secondary", ]
dt_both = dt [which(dt$V3 == "both_fp" & dt$V1 == "primary"), ]

total_reads_p = length(unique(dt_p$V4))
grob_read_count_p = grobTree(textGrob(paste0("Total Reads: ", total_reads_p),
                                    x = 0.85, y = 0.95, hjust = 0.5))
total_reads_s = length(unique(dt_s$V4))
grob_read_count_s = grobTree(textGrob(paste0("Total Reads: ", total_reads_s),
                                    x = 0.85, y = 0.95, hjust = 0.5))
total_reads_u = length(unique(dt$V4)) 
grob_read_count_u = grobTree(textGrob(paste0("Total Reads: ", total_reads_u),
                                    x = 0.85, y = 0.95, hjust = 0.5))

total_reads_both = length(unique(dt_both$V4))
grob_read_count_both = grobTree(textGrob(paste0("Total Reads: ", total_reads_both),
                                    x = 0.85, y = 0.95, hjust = 0.5))

p_title = paste0(args[2], " (primary peak)") 
s_title = paste0(args[2], " (secondary peak at ", args[3], ")") 
both_title = paste0(args[2], " - footprints covering primary & secondary") 

x_tic_breaks = seq(0,300,  5)
x_tic_breaks_label = seq(0,300,  25)
plt_p = ggplot(dt_p, aes(x = V2)) + 
        scale_x_continuous(name = "Footprint length [bp]", 
                         breaks = x_tic_breaks_label, labels = x_tic_breaks_label,
                         limits = c(0,300)) + 
        geom_histogram(position = "identity",  breaks = x_tic_breaks) + 
        geom_rangeframe() + theme_few() + ggtitle(p_title) + 
        xlab("Footprint length [bp]") + theme (plot.title = element_text(hjust = 0.5)) + annotation_custom(grob_read_count_p) 
        #scale_x_continuous(name = "Footprint length [bp]", breaks = x_tic_breaks, limits = c(0, max(dt_p$V2)))
plt_s = ggplot(dt_s, aes(x = V2)) +  
        scale_x_continuous(name = "Footprint length [bp]", 
                         breaks = x_tic_breaks_label, labels = x_tic_breaks_label,
                         limits = c(0,300)) + 
        geom_histogram(position = "identity", breaks = x_tic_breaks) + 
        geom_rangeframe() + theme_few() + ggtitle(s_title) + 
        xlab("Footprint length [bp]") + theme (plot.title = element_text(hjust = 0.5)) + annotation_custom(grob_read_count_s)

plt_both = ggplot(dt, aes(x = V2, fill = V1)) + 
        geom_histogram(position = "identity", bins = 30, alpha = 0.4) + 
        geom_rangeframe() + theme_few() + ggtitle(paste0(args[2] ," - layover")) + 
        theme(legend.position = c(0.8, 0.75)) + 
        xlab("Footprint length [bp]") + theme (plot.title = element_text(hjust = 0.5)) + annotation_custom(grob_read_count_u) + 
        labs(fill = "Peak type") 
plt_p_and_s = ggplot(dt_both, aes(x = V2)) +  
        geom_histogram(position = "identity", bins = 30) + 
        geom_rangeframe() + theme_few() + ggtitle(both_title) + 
        xlab("Footprint length [bp]") + theme (plot.title = element_text(hjust = 0.5)) + annotation_custom(grob_read_count_both)


plt = cowplot::plot_grid(plt_p, plt_s, plt_both, plt_p_and_s, nrow = 2, ncol = 2)
Cairo::CairoPNG(args[4], width = 10, height = 8, res = 150, units = "in")
print(plt)
dev.off() 

pdf(args[5], width = 10, height = 8)
print (plt)
dev.off()

