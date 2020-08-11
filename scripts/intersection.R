suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("IRanges"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("yaml"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("stringi"))
suppressPackageStartupMessages(library("R.utils"))
source("scripts/utils.R")
options(traceback = T)
#options(error=function() { traceback(2); if(!interactive())

args = commandArgs(trailingOnly = T)

target_bed = args[1]
methylation_vec_bedgz = args[2]
chunk_size = as.integer(args[3])
out_fp = file(args[4], "w")

target_bed_dt = read.csv(target_bed, header = F, sep = "\t", 
                         stringsAsFactors = F, comment.char = "#")
names(target_bed_dt) = bed_file_header
target_bed_dt$`chromStop` = target_bed_dt$`chromStop` + 1 

#target_bed_gt = GenomicRanges::GNCList(target_bed_gr)
print (head(target_bed_dt))
methylation_vec_con = NULL

if (is_gz(methylation_vec_bedgz)){
    methylation_vec_con = file(methylation_vec_bedgz, "rt")
}else{
    methylation_vec_con = file(methylation_vec_bedgz, "r")
}
while_loop_counter = 0
valid_fragments_count = 0
while (TRUE){
    while_loop_counter = while_loop_counter + 1
    methylation_chunk = read_from_connection(con = methylation_vec_con, 
                                chunk_size = chunk_size)
    if (!is.null(methylation_chunk)){
        names(methylation_chunk) = methylated_bed_header
        # creat interval tree for methylation chunk
        methylation_chunk_gr = GRanges(IRanges(start = methylation_chunk$`m_chromStart`, 
                                               end = methylation_chunk$`m_chromStop`,),
                                       seqnames = methylation_chunk$`m_chrom`)
        methylation_chunk_gt = GenomicRanges::GNCList(methylation_chunk_gr)
       
        overlap.idx = fragments_covering_enhancer_peak (
                         df = target_bed_dt, it = methylation_chunk_gt)
        target_bed_sub_dt = target_bed_dt[overlap.idx@from, ]
        methylation_chunk = methylation_chunk[overlap.idx@to, ]
        if (nrow(target_bed_sub_dt)!=0){
            intersect_df = cbind(target_bed_sub_dt, methylation_chunk)
            #print (head(intersect_df))
            # #################### Not Required ####################
            # reverse the string to account for peak strandedness
            #strand_aware_vec = unlist(lapply(seq(nrow(intersect_df)), 
            #  function (x){
            #      if (identical(intersect_df[x, "strand"], "-")){
            #         return(stringi::stri_reverse(intersect_df[x, "m_vec"] ))
            #      }else{
            #         return(intersect_df[x, "m_vec"])
            #      }
            #  }))
            #intersect_df["m_vec"] = strand_aware_vec 
            #########################################################
            write.table(intersect_df, out_fp, col.names = F, row.names = F,
                        quote = F, sep = "\t", append = T)

        }else{
            next
        }
        
    }else{
        cat ("Finished reading from file", 
             c(methylation_vec_bedgz), 
             "closing connection!",  "\n")
        #close(cfdna_con)
        break
    }
}

