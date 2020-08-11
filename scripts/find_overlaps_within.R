suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("IRanges"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("yaml"))
suppressPackageStartupMessages(library("stringr")) 
source("scripts/utils.R")
options(show.error.locations = T)
args = commandArgs(trailingOnly = T)
mapped_read_bed = args[2]
out_con = file(args[3], "w")
peak.data = read.table(args[1], sep = "\t", header = F,
           stringsAsFactors = F, comment.char = "#")
names(peak.data) = vplot_peak_pair_header
#peak.gr = GRanges(IRanges(start = peak.data$chromStart, end = peak.data$chromStop),
#           seqnames = peak.data$chrom)
#peak.git = GenomicRanges::GNCList(peak.gr)

mapped_con = file(mapped_read_bed, "r") 

while_loop_counter = 0
mapped_line_count = as.double(0)
valid_fragments_count = 0
while (TRUE){
    while_loop_counter = while_loop_counter + 1
    mapped_chunk = read_from_connection(con = mapped_con, 
        chunk_size = 10000)
    if (!is.null(mapped_chunk)){
        names(mapped_chunk) = mapped_read_header
        #print (head(mapped_chunk))
        #print (head(peak.data))
        mapped_chunk_gr = GRanges(IRanges(start = mapped_chunk$chromStart,
                                          end = mapped_chunk$chromStop),
                                  seqnames = mapped_chunk$chrom)
        mapped_chunk_git = GenomicRanges::GNCList(mapped_chunk_gr)
        hits = fragments_covering_enhancer_range(df = peak.data, it = mapped_chunk_git)
        #print (c(length(hits@from), length(hits@to)))
        peak_intersect_df = peak.data[hits@from, ]
        
        mapped_reads_intersect_df = mapped_chunk[hits@to, ]
        #row.names(peak_intersect_df) = seq(nrow(peak_intersect_df))
        #row.names(mapped_reads_intersect_df) = seq(nrow(mapped_reads_intersect_df))
        #print (c(nrow(peak_intersect_df), nrow(mapped_reads_intersect_df)))
        if (nrow(mapped_reads_intersect_df) !=0){ 
            intersect_df = cbind (peak_intersect_df, mapped_reads_intersect_df)
            write.table(intersect_df, out_con, col.names = F, row.names = F, 
                        quote = F, sep = "\t", append = T)
            #print (while_loop_counter*1000)
        }else{
            next
        }
    }else{
        cat ("Finished reading from file", 
             c(mapped_read_bed), 
             "closing connection!",  "\n")
        #close(cfdna_co)
        break
    }
}

close(out_con)
