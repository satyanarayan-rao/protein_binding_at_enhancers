library(stringr)
library(ggplot2)
library(ggthemes)
library(zoo)

# chr2L@11806479@11806979@chr2L:11806729-11806729^20
args = commandArgs(trailingOnly = T)
# args[1]: mnase csv gz file
# args[2]: site of interest 
# args[3]: site wise eom plot
# args[4]: site wise erom rmean plot
# args[5]: site wise max peak location (have to find two peaks)
# args[6]: rolling mean window size
input_dt = read.table(args[1], stringsAsFactors = F, sep = "", header = F)
site_of_interest = str_replace(args[2], "\\^", "\\\\^") 
rolling_window_size = as.integer (args[6])
# get the specific site
site_dt = as.data.frame(t(input_dt[grepl(site_of_interest, input_dt$V1),, drop = F]),
                  stringsAsFactors = F)
# remove the first row
site_dt = site_dt[2:dim(site_dt)[1], , drop = F]
names(site_dt) = "site_of_interest"
site_dt[site_of_interest] = as.numeric(site_dt[, "site_of_interest"])
site_dt["x_tics"] = seq(dim(site_dt)[1])
site_dt$x_tics = site_dt$x_tics - as.integer(dim(site_dt)[1]/2)
site_dt$rmean = rollmean(site_dt[site_of_interest], 
                  k = rolling_window_size, fill = NA) 

# find the max (indices of right and left flank) and the second max
###

to_find_max = site_dt[, "site_of_interest"]
total_rows = dim(site_dt)[1]

first_max_loc = which.max(to_find_max)
first_max_val = to_find_max[first_max_loc]
first_max_right_flank = NULL
for (i in seq(first_max_loc + 1, total_rows)){
    if (to_find_max[i] != first_max_val){
        first_max_right_flank = i - 1
        break
    }else{
        next
    }
}
first_max_left_flank = NULL
for (i in seq(first_max_loc - 1 , 1)){
    if (to_find_max[i] != first_max_val){
        first_max_left_flank = i + 1
        break
    }else{
        next
    }
}

# copy to_find_max to another vector and replace first_max val to -1 and then find max
print (c (first_max_left_flank, first_max_right_flank))
to_find_second_max = to_find_max
to_find_second_max[first_max_left_flank: first_max_right_flank] = -1 
second_max_loc =  which.max (to_find_second_max)
second_max_val = to_find_second_max[second_max_loc] 
print (c(second_max_val, second_max_loc) )
second_max_right_flank = NULL 
for (i in seq(second_max_loc + 1, total_rows)){
    print(to_find_second_max[i])
    if (to_find_second_max[i] != second_max_val){
        second_max_right_flank = i - 1
        break
    }else{
        next
    }
}
print ('here')
second_max_left_flank = NULL
for (i in seq(second_max_loc - 1 , 1)){
    if (to_find_second_max[i] != second_max_val){
        second_max_left_flank = i + 1
        break
    }else{
        next
    }
}
first_max_left_loc_eom = site_dt[first_max_left_flank, "x_tics"]
first_max_right_loc_eom = site_dt[first_max_right_flank, "x_tics"]
second_max_left_loc_eom = site_dt[second_max_left_flank, "x_tics"]
second_max_right_loc_eom = site_dt[second_max_right_flank, "x_tics"]

###
# find the max (indices of right and left flank) and the second max: for rolling mean
###

to_find_max = site_dt[, "rmean"]
total_rows = dim(site_dt)[1]

first_max_loc = which.max(to_find_max)
first_max_val = to_find_max[first_max_loc]
first_max_right_flank = NULL
for (i in seq(first_max_loc + 1, total_rows)){
    if (to_find_max[i] != first_max_val){
        first_max_right_flank = i - 1
        break
    }else{
        next
    }
}
first_max_left_flank = NULL
for (i in seq(first_max_loc - 1 , 1)){
    if (to_find_max[i] != first_max_val){
        first_max_left_flank = i + 1
        break
    }else{
        next
    }
}

# copy to_find_max to another vector and replace first_max val to -1 and then find max
to_find_second_max = to_find_max
to_find_second_max[first_max_left_flank: first_max_right_flank] = -1 
second_max_loc =  which.max (to_find_second_max)
second_max_val = to_find_second_max[second_max_loc] 
second_max_right_flank = NULL 
for (i in seq(second_max_loc + 1, total_rows)){
    if (to_find_second_max[i] != second_max_val){
        second_max_right_flank = i - 1
        break
    }else{
        next
    }
}
second_max_left_flank = NULL
for (i in seq(second_max_loc - 1 , 1)){
    if (to_find_second_max[i] != second_max_val){
        second_max_left_flank = i + 1
        break
    }else{
        next
    }
}
first_max_left_loc_rmean = site_dt[first_max_left_flank, "x_tics"]
first_max_right_loc_rmean = site_dt[first_max_right_flank, "x_tics"]
second_max_left_loc_rmean = site_dt[second_max_left_flank, "x_tics"]
second_max_right_loc_rmean = site_dt[second_max_right_flank, "x_tics"]


###



pdf(args[3], width = 6, height = 6)
kk = ggplot (data = site_dt, aes(x = x_tics, y = site_of_interest)) + 
         geom_line() + ggtitle(site_of_interest) + theme_few() + 
         geom_rangeframe() + theme(plot.title = element_text(hjust = 0.5)) + 
         xlab("Distance from peak center [bp]") + ylab ("MNase (e.o.m)") + 
         geom_vline(xintercept = c(first_max_left_loc_eom, first_max_right_loc_eom,
                                   second_max_left_loc_eom, second_max_right_loc_eom),
                    linetype = 2)


print(kk)
dev.off()

pdf(args[4], width = 6, height = 6)
kk = ggplot (data = site_dt, aes(x = x_tics, y = rmean)) + 
         geom_line() + ggtitle(site_of_interest) + theme_few() + 
         geom_rangeframe() + theme(plot.title = element_text(hjust = 0.5)) + 
         xlab("Distance from peak center [bp]") + ylab ("MNase (e.o.m-rolling mean)") + 
         geom_vline(xintercept = c(first_max_left_loc_rmean, first_max_right_loc_rmean,
                                   second_max_left_loc_rmean, second_max_right_loc_rmean),
                    linetype = 2)

print(kk)
dev.off()

# save max locations in file and also minimum difference between peaks flanks
# create a data frame of locations and differences
location_df = data.frame(row.names = c("eom", "roll_mean"),
                first_max_left_loc = c(first_max_left_loc_eom, 
                                       first_max_left_loc_rmean),
                first_max_right_loc = c(first_max_right_loc_eom, 
                                        first_max_right_loc_rmean),
                second_max_left_loc = c(second_max_left_loc_eom,
                                        second_max_left_loc_rmean),
                second_max_right_loc = c(second_max_right_loc_eom,
                                         second_max_right_loc_rmean))
write.table(location_df, file = args[5], row.names = T, col.names = T, quote = F)
