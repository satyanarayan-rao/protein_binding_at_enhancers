library(ggplot2)
library(ggthemes)
library(Cairo)
library(grid)
args = commandArgs(trailingOnly = T)

# ------------------- #



# ------------------- #

figure2_theme <- function (){
    theme (axis.text.y =element_text(vjust =1))+
    theme(plot.title=element_text( size=20 )) +
    theme(axis.title.x = element_text(colour = "black", size = 20),
          axis.title.y = element_text(colour = "black", size = 20)) +
    theme(axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15)) +
    theme(legend.title= element_text(size = 20),
          legend.text = element_text(size = 20)) +
    theme(axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5))
}
 
cut_str = paste0("cut -f", args[3], " ", args[1])
per_a = read.table(pipe(cut_str)) 
cut_str = paste0("cut -f", args[3], " ", args[2])
per_b = read.table(pipe(cut_str)) 

bw = as.integer(args[4])
per_a$V1 = (per_a$V1)*100
per_b$V1 = (per_b$V1)*100

fixed_breaks = seq(0, 100, bw)
hist_a = hist(per_a$V1, breaks = fixed_breaks, plot = F) 
hist_b = hist(per_b$V1, breaks = fixed_breaks, plot = F) 

dens_df_a = data.frame(dens = hist_a$density, x_tics = hist_a$mids, label = "Open") 
dens_df_b = data.frame(dens = hist_b$density, x_tics = hist_b$mids, label = "Closed")
nrow_a = nrow (per_a)
nrow_b = nrow (per_b)
color_open_and_close = c("Open" = "#386cb0", "Closed" = "#f0027f")

labels_vec = c (paste0("Open (#reads = ", format(nrow_a, big.mark = ",", scientific = F), ")"),
                paste0("Closed (#reads = ", format(nrow_b, big.mark=",", scientifc=F), ")"))
#nrow_a_text = grobTree(textGrob(paste0("Mapped reads to open (n = ", nrow_a, ")"), x = 0.5, y = 0.85, hjust = 0.5))
#nrow_b_text = grobTree(textGrob(paste0("Mapped reads to closed (n = ", nrow_b, ")"), x = 0.5, y = 0.8, hjust = 0.5))  

plot_df = rbind(dens_df_a, dens_df_b)

pdf(args[5], width = 6, height = 4)
plt = ggplot(plot_df, aes (x = x_tics, y = dens, color = label)) +
      geom_line() + geom_rangeframe() + theme_few() + 
      xlab("Percentage methylation") + ylab("Frequency [A.U.]") + 
     # annotation_custom(nrow_a_text) + annotation_custom(nrow_b_text) + 
      scale_color_manual(values = color_open_and_close, 
                         breaks = c("Open", "Closed"),
                         labels = c(labels_vec)) + 
      theme(legend.position = c(0.7, 0.7), legend.text.align = 1) + 
      labs(color = "") + 
      scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + 
      figure2_theme() 

print(plt)
dev.off()

Cairo::CairoPNG(args[6], width = 6, height = 4, units = "in", res = 150)
print(plt)
dev.off()

############# do smoothing using spline #################

sp_hist_a = smooth.spline(x = hist_a$mids , y = hist_a$density, df = 10)
sp_hist_b = smooth.spline(x = hist_b$mids , y = hist_b$density, df = 10)

pred_a = predict(sp_hist_a, x = hist_a$mids)
pred_y_a = pred_a$y
pred_b = predict(sp_hist_b, x = hist_b$mids)
pred_y_b = pred_b$y

dens_df_a = data.frame(dens = pred_y_a, x_tics = hist_a$mids, label = "Open") 
dens_df_b = data.frame(dens = pred_y_b, x_tics = hist_b$mids, label = "Closed")

plot_df = rbind(dens_df_a, dens_df_b)

pdf(args[7], width = 6, height = 4)
plt = ggplot(plot_df, aes (x = x_tics, y = dens, color = label)) +
      geom_line() + geom_rangeframe() + theme_few() + 
      xlab("Percentage methylation") + ylab("Frequency [A.U.]") + 
      #annotation_custom(nrow_a_text) + annotation_custom(nrow_b_text) +
      scale_color_manual(values = color_open_and_close, 
                         breaks = c("Open", "Closed"),
                         labels = c(labels_vec)) + 
      theme(legend.position = c(0.7, 0.7), legend.text.align = 1) + 
      labs(color = "")



print(plt)
dev.off()

Cairo::CairoPNG(args[8], width = 6, height = 4, units = "in", res = 150)
print(plt)
dev.off()

#########################################################



################ Do a emperical CDF just by doing cumsum #########  

cumsum_df_a = data.frame(x_tics = hist_a$mids, cdf = cumsum(pred_y_a)/max(cumsum(pred_y_a)), label = "Open")
cumsum_df_b = data.frame(x_tics = hist_b$mids, cdf = cumsum(pred_y_b)/max(cumsum(pred_y_b)),
 label = "Closed")

plot_df = rbind(cumsum_df_a, cumsum_df_b) 

pdf(args[9], width = 6, height = 4)
plt = ggplot(plot_df, aes (x = x_tics, y = cdf, color = label)) +
      geom_line() + geom_rangeframe() + theme_few() + 
      xlab("Percentage methylation") + ylab("eCDF") + 
      #annotation_custom(nrow_a_text) + annotation_custom(nrow_b_text) +
      scale_color_manual(values = color_open_and_close, 
                         breaks = c("Open", "Closed"),
                         labels = c(labels_vec)) + 
      theme(legend.position = c(0.75, 0.25), legend.text.align = 1) + 
      labs(color = "") + geom_vline(xintercept = 50, lty = 2)



print(plt)
dev.off()

Cairo::CairoPNG(args[10], width = 6, height = 4, units = "in", res = 150)
print(plt)
dev.off()

########################### normalize count based eCDF #############
dens_df_a = data.frame(dens = (hist_a$counts)/nrow_a, x_tics = hist_a$mids, label = "Open") 
dens_df_b = data.frame(dens = (hist_b$counts)/nrow_b, x_tics = hist_b$mids, label = "Closed") 

plot_df  = rbind(dens_df_a, dens_df_b)
print(plot_df)

pdf(args[11], width = 6, height = 4)
plt = ggplot(plot_df, aes (x = x_tics, y = dens, color = label)) +
      geom_line() + geom_rangeframe() + theme_few() + 
      xlab("Percentage methylation") + ylab("Frequency [A.U.]") + 
     # annotation_custom(nrow_a_text) + annotation_custom(nrow_b_text) + 
      scale_color_manual(values = color_open_and_close, 
                         breaks = c("Open", "Closed"),
                         labels = c(labels_vec)) + 
      theme(legend.position = c(0.7, 0.7), legend.text.align = 1) + 
      labs(color = "")

print(plt)
dev.off()

Cairo::CairoPNG(args[12], width = 6, height = 4, units = "in", res = 150)
print(plt)
dev.off()


