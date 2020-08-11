processFile = function(filepath, chunk_size = 1000) {
  con = file(filepath, "r")
  while ( TRUE ) {
    line = readLines(con, n = chunk_size)
    if ( length(line) == 0 ) {
      break
    }
    dat = read.csv (text=line, sep="\t", header=FALSE)
    print(dat)
  }

  close(con)
}

read_from_connection = function(con=NULL, 
                                chunk_size=1, 
                                sep="\t"){
    lines = readLines(con, n = chunk_size) 
    #cat (paste0(length(lines), "\n"))
    if(length(lines) == 0){
	return (NULL) # cat("Finished reading from connection!\n")
    }else{
        df = read.csv(text=lines, sep=sep, header=FALSE, stringsAsFactors = F)
        return(df)
    }
}

create_tfbs_hash = function (df = NULL, seed_count = NULL, label = label){
    # df: tfbs df
    # seed_count: usually a vector of zeros 
    # label: label for seed count
    keys = do.call(paste, c(df[, seq(4)], sep="^"))
    tfbs_dict = list() 
    label_and_count = setNames(as.list(seed_count), label)
    for (k in keys){
        tfbs_dict[[k]] = label_and_count
    }
    return (tfbs_dict)
}

fragments_covering_enhancer_peak = function(df=NULL, it=NULL){
    # df: data frame; obtained from pairs file 
    # it: interval tree generated from bed file
    # reads: a list containing parameters to process reads; 
    #        for example, length constrain, fill_in reads. 
    #        see configuration yaml for details 
    # returns: findOverlaps object 
    #mid_point = as.integer((df$chromStop - df$chromStart)/2 + 0.5) 
    # create an GRanges object from df
    df.gr = GRanges (IRanges(start = df$`chromStart`, 
                             end = df$`chromStop`),
                     seqnames = df$`chrom`)
    
    # find overlap with TFBS
    #overlapped_reads = findOverlaps(df.gr, it, type="within", select = "all")
    overlapped_reads = findOverlaps(df.gr, it, select = "all")
    
    return (overlapped_reads)
    
}

get_other_file_name = function (f, prefix="count", sep="_"){
    
    splitted = unlist(strsplit(f, "/"))
    splitted[length(splitted)] = paste (
                  prefix, 
                  splitted[length(splitted)], 
                  sep=sep)
    new_file_path = paste (splitted, collapse="/")
    return (new_file_path)
     

}

pairs_file_header= c ("chrom", 
                      "sam_flag",
                      "chromStart",
                      "chromStop",
                      "arbit_score",
                      "fragment_length")

bed_file_header= c ("chrom", 
                    "chromStart",
                    "chromStop",
                    "name",
                    "score",
                    "strand")
methylated_bed_header = c ("m_chrom", 
                    "m_chromStart",
                    "m_chromStop",
                    "m_name",
                    "m_score",
                    "m_vec")
# example line
# 
# chr3R	19747738	19748134	peak_3545	chr3R:19748119-19748119^3	-	0	366
vplot_peak_pair_header = c(
                    "chrom", 
                    "chromStart",
                    "chromStop",
                    "peak_id",
                    "chr_loc",
                    "strand",
                    "primary_peak",
                    "secondary_peak"
)
# mapped methylated reads
# chr2L	11806460	11806747	SRR3133326.952549_952549/1_overlapping`99~147	.
mapped_read_header = c(
                    "chrom", 
                    "chromStart",
                    "chromStop",
                    "read_id",
                    "name"
)

############## find overalps where enahncer range completely fall within read ############# 
fragments_covering_enhancer_range = function(df=NULL, it=NULL){
    df.gr = GRanges (IRanges(start = df$`chromStart`, 
                             end = df$`chromStop`),
                     seqnames = df$`chrom`)
    
    overlapped_reads = findOverlaps(df.gr, it, select = "all", type = "within")
    
    return (overlapped_reads)
}





is_gz = function (filename){
    details = summary(file(filename))
    if (details$class == "gzfile") {
        return (TRUE)
    } else {
        return (FALSE)
    }
}
# Not using `init_enrichment_dt` anymore - seems like populating a dataframe and doing operations
# (like incrementing cell values is too slow)

#init_enrichment_dt = function(roi_df=NULL){
#    # roi_df: a bed dataframe. 
#    # return: a data.table with chrom loci as rownames and relative offsets("as
#    # character") as column names 
#    roi_max_length = max (roi_df$chromStop - roi_df$chromStart + 1 )
#    #print (roi_max_length)
#    reads_enrichment_dt = data.frame( # its important to initialize a data.table, not data.frame 
#                                     matrix (0L, 
#                                     nrow = nrow(roi_df),
#                                     ncol = roi_max_length))
#    rownames_dt = paste(
#                        roi_df$chrom, 
#                        roi_df$chromStart, 
#                        roi_df$chromStop, sep ="_")
#    rownames(reads_enrichment_dt) = rownames_dt
#    rowname_to_index_map = data.frame ("idx" = seq(1, nrow(reads_enrichment_dt)),
#                                       row.names = rownames_dt)
#    if (roi_max_length %% 2 == 1){ # if the roi_length is odd
#        relative_pos = as.character(seq(0 - (roi_max_length - 1)/2, 
#                                        0 + (roi_max_length - 1)/2))
#        names(reads_enrichment_dt) = relative_pos
#       
#    } else{ 
#        relative_pos = as.character(seq(0 - roi_max_length/2, 
#                                        0 + roi_max_length/2 - 1))
#        names(reads_enrichment_dt) = relative_pos
#    } 
#    colname_to_index_map = data.frame (idx = seq(1, ncol(reads_enrichment_dt)),
#                                       row.names = names(reads_enrichment_dt))
#
#    return (list ("reads_enrichment_dt" = reads_enrichment_dt,
#                  "rowname_to_index_map" = rowname_to_index_map,
#                  "colname_to_index_map" = colname_to_index_map)) 
#} 

# Not using `fill_in_reads` function anymore - experienced that its too slow 

#fill_in_reads = function (reads_enrichment_dt = NULL, 
#                          roi_df = NULL,
#                          reads_df = NULL,
#                          fill_in_region = NULL
#                          #rowname_to_index_map = NULL,
#                          #colname_to_index_map = NULLi
#                          ){
#    # reads_enrichment_dt: a data.table that will be updated inplace
#    # roi_df: a datframe where paired reads have overlapped 
#    # reads_df: a paired read dataframe that mapped to roi
#    # fill_in_region: number of bases to be enriched (1 added) upstream and
#    # downstream 
#    # returns: the modified reads_enrichment_dt 
#    row_indices = paste (roi_df$chrom, roi_df$chromStart, roi_df$chromStop, sep="_")
#    left_flank_relative_pos = as.integer (names(reads_enrichment_dt)[1]) # the first column name indicates of left upstream
#    right_flank_relative_pos = as.integer (names(reads_enrichment_dt)[-1]) # the first column name indicates of left upstream
#    mid_point = left_flank_relative_pos + 
#                (reads_df$chromStart - roi_df$chromStart) + 
#                as.integer(reads_df$fragment_length/2 + 0.5)    
#    # create a data frame of row_indices, mid_point, and fragment lengths 
#    indices_to_update = data.frame ("row_indices" = row_indices, 
#                                    "mid_point" = mid_point, 
#                                    "fragment_length" = reads_df$fragment_length,
#                                    stringsAsFactors = FALSE)  
#    for (i in seq(1, nrow(indices_to_update))){
#        row_indices = indices_to_update[i, "row_indices"]
#        mid_point = indices_to_update[i, "mid_point"]
#        fragment_length = indices_to_update[i, "fragment_length"]
#        left_flank = min(fill_in_region$upstream, fragment_length/2) 
#        right_flank = min(fill_in_region$downstream, fragment_length/2)
#        for (j in seq(max (left_flank_relative_pos, mid_point - left_flank),
#                      min (right_flank_relative_pos, mid_point + right_flank))){
#            #reads_enrichment_dt[as.integer(rowname_to_index_map[row_indices, "idx"]), as.character(j)] = reads_enrichment_dt[as.integer(rowname_to_index_map[row_indices, "idx"]), as.character(j)] + 1 
#            
#            reads_enrichment_dt[row_indices, as.character(j)] = reads_enrichment_dt[row_indices, as.character(j)] + 1 
#            #reads_enrichment_dt[2, "-1005"] = reads_enrichment_dt[2, "-1005"] + 1 
#        }
#        
#    }
#    return (reads_enrichment_dt)
#}


