################################## co-binding analysis ##################################### 

rule extend_primary_peak_for_cobinding_analysis: 
    input:
        primary_peak_data = "peak_footprints/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_sam_flag_{sam_flag}.tsv",
        footprint_capped_pkl = "footprint_min_size_selection/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_tf_footprint_min_length_10_wobble_gap_1.pkl"
    params:
        strand = lambda wildcards: get_strand_info(wildcards)
        
    output:
        cobinding_footprint_vec = "cobinding_footprints/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_lex_{lextend}_rex_{rextend}_site_{site}_sam_flag_{sam_flag}.tsv",
        cobinding_footprint_vec_verbose = "cobinding_footprints/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_lex_{lextend}_rex_{rextend}_site_{site}_sam_flag_{sam_flag}_verbose.tsv"
    shell:  
        "python scripts/extend_footprint.py {input.primary_peak_data}"
        " {input.footprint_capped_pkl} {wildcards.lflank} {wildcards.rflank}"
        " {wildcards.lextend} {wildcards.rextend}"
        " {output.cobinding_footprint_vec}"
        " {params.strand}"
        " {output.cobinding_footprint_vec_verbose}" 

rule select_reads_for_cobinding_analysis:
    input:
        cobinding_footprint_vec = "cobinding_footprints/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_lex_{lextend}_rex_{rextend}_site_{site}_sam_flag_{sam_flag}.tsv",
        #primary_secondary_read_map_pkl = "vplot_lists/intersected_reads.pkl"
        methylated_read_start_and_end_dict = "footprints/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_tf_footprint.read_dict.pkl"
    params:
        strand = lambda wildcards: get_strand_info(wildcards)
    output:
        filtered_cobinding_vec = "filtered_cobinding_reads/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_sec_peak_at_{s_peak}_lex_{lextend}_rex_{rextend}_site_{site}_sam_flag_{sam_flag}.tsv",
        
    shell:
        "python scripts/select_reads_for_cobinding_analysis.py" 
        " {input.cobinding_footprint_vec} {input.methylated_read_start_and_end_dict}"
        " {wildcards.site} {wildcards.s_peak} {params.strand}"
        " {wildcards.sam_flag} {output.filtered_cobinding_vec}"
        " {wildcards.lflank} {wildcards.rflank}"

rule fplen_dist_at_primary_and_sec_peaks: 
    input:
        filtered_cobinding_vec = "filtered_cobinding_reads/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_sec_peak_at_{s_peak}_lex_{lextend}_rex_{rextend}_site_{site}_sam_flag_{sam_flag}.tsv",
        occluded_pkl = "occluded_edges_on_methylation_vec/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_min_flen_10_wobble_gap_1_occluded.pkl", 
     
    output:
        hist_tsv = "fplen_filtered_cobinding_reads/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_sec_peak_at_{s_peak}_lex_{lextend}_rex_{rextend}_site_{site}_sam_flag_{sam_flag}.fplen.tsv"

    shell:
        "python scripts/get_fplen_at_primary_and_secondary.py"
        " {input.filtered_cobinding_vec} {wildcards.lflank} {wildcards.rflank}"
        " {wildcards.s_peak} {wildcards.lextend} {output.hist_tsv}"
        " {input.occluded_pkl}"

rule plot_fplen_primary_and_secondary_peaks: 
    input:
        hist_tsv = "fplen_filtered_cobinding_reads/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_sec_peak_at_{s_peak}_lex_{lextend}_rex_{rextend}_site_{site}_sam_flag_{sam_flag}.fplen.tsv"
    params:
    output:
        fplen_png = "plots/fplen_dist_primary_and_secondary/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_sec_peak_at_{s_peak}_lex_{lextend}_rex_{rextend}_site_{site}_sam_flag_{sam_flag}.fplen.png",
        fplen_pdf = "plots/fplen_dist_primary_and_secondary/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_sec_peak_at_{s_peak}_lex_{lextend}_rex_{rextend}_site_{site}_sam_flag_{sam_flag}.fplen.pdf"
    shell:
        "Rscript scripts/plot_fplen_primary_and_secondary.R"
        " {input.hist_tsv} {wildcards.site} {wildcards.s_peak}"
        " {output.fplen_png} {output.fplen_pdf}"

rule assign_cobinding_modes:
    input:
        filtered_cobinding_vec = "filtered_cobinding_reads/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_sec_peak_at_{s_peak}_lex_{lextend}_rex_{rextend}_site_{site}_sam_flag_{sam_flag}.tsv",
        occluded_pkl = "occluded_edges_on_methylation_vec/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_min_flen_10_wobble_gap_1_occluded.pkl", 
    params:
    output:
        annotated_vec = "annotated_cobinding_status/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_sec_peak_at_{s_peak}_lex_{lextend}_rex_{rextend}_site_{site}_sam_flag_{sam_flag}_cobinding.tsv", 
        annotated_vec_150 = "annotated_cobinding_status/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_sec_peak_at_{s_peak}_lex_{lextend}_rex_{rextend}_site_{site}_sam_flag_{sam_flag}_150_cobinding.tsv", 
    shell: 
        "python scripts/assign_cobinding.py {input.filtered_cobinding_vec}"
        " {wildcards.lflank} {wildcards.rflank} {wildcards.s_peak}"
        " {wildcards.lextend} {output.annotated_vec} {input.occluded_pkl}"
        " {output.annotated_vec_150}" 

def get_files_for_excess_ratio_cobinding(wildcards): 
    flist = [] 
    peak_fp = open(config[wildcards.setting])
    for line in peak_fp:
        line_items = line.strip().split()
        peak_name = line_items[0]
        cl_id = line_items[1]  
        s_peak = line_items[2]
        fname = "annotated_cobinding_status/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_sec_peak_at_{s_peak}_lex_{lextend}_rex_{rextend}_site_{site}_sam_flag_{sam_flag}_cobinding.tsv".format(
           sample = wildcards.sample, 
           bed = wildcards.bed,
           cl_id = cl_id,
           lflank = wildcards.lflank,
           rflank = wildcards.rflank, 
           site = peak_name, 
           sam_flag = wildcards.sam_flag,
           s_peak = s_peak,
           lextend = wildcards.lextend, 
           rextend = wildcards.rextend) 
        flist.append(fname)
    return flist    
rule excess_ratio_for_cobinding:
    input:
        all_files = lambda wildcards: get_files_for_excess_ratio_cobinding (wildcards)
    params: 
    output:
        excess_ratio_flist_tsv = "flist_excess_ratio_cobinding/{sample}_to_{bed}_excess_ratio_for_{setting}_with_lf_{lflank}_rf_{rflank}_lex_{lextend}_rex_{rextend}_sam_flag_{sam_flag}_excess_ratio.flist.tsv"
    run:
        fp = open (output.excess_ratio_flist_tsv, "w")
        for f in input.all_files:
            fp.write(f + "\n")
        fp.close()

rule get_excess_ratio: 
    input:
        excess_ratio_flist_tsv = "flist_excess_ratio_cobinding/{sample}_to_{bed}_excess_ratio_for_{setting}_with_lf_{lflank}_rf_{rflank}_lex_{lextend}_rex_{rextend}_sam_flag_{sam_flag}_excess_ratio.flist.tsv"
    params:
    output:
        excess_ratio_tsv = "excess_ratio_cobinding/{sample}_to_{bed}_excess_ratio_for_{setting}_with_lf_{lflank}_rf_{rflank}_lex_{lextend}_rex_{rextend}_sam_flag_{sam_flag}_excess_ratio.tsv"
    shell:
        "python scripts/excess_ratio_cobinding.py {input.excess_ratio_flist_tsv}"
        " {output.excess_ratio_tsv}"
#rule percentage_occupancy_cobinding:
     
