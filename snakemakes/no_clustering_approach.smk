rule select_for_site:
    input:
        methylation_matrix_for_cl =  "flank_footprint_matrix/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_methylation_matrix.cluster_{cl_id}.tsv"
    params:
        string_for_site = lambda wildcards: config["site_annotation"][wildcards.cl_id][wildcards.site] + "\^" + wildcards.cl_id
    output:
        peak_data = "peak_footprints/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_sam_flag_{sam_flag}.tsv", 
    shell: 
        "grep \"{params.string_for_site}\" {input.methylation_matrix_for_cl} | grep -w {wildcards.sam_flag} > {output.peak_data}"


rule get_real_footprint_lengths_and_orange_on_read: 
    input:
        peak_data = "peak_footprints/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_sam_flag_{sam_flag}.tsv", 
        methylation_matrix_real_footprint_len_pkl =\
              "flank_footprint_matrix/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_methylation_matrix_flen.pkl",
        footprint_capped_pkl = "footprint_min_size_selection/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_tf_footprint_min_length_10_wobble_gap_1.pkl"
    params:
        neighbor_left_from_peak_center = 30,
        neighbor_right_from_peak_center = 30,
        strand = lambda wildcards: get_strand_info(wildcards)
        
    output:
        read_per_orange_and_length = "peak_per_orange_vs_footprint_length/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_po_vs_footplen_for_site_{site}_sam_flag_{sam_flag}.tsv", # po: percentage orange; footplen: footprint length
        read_per_methylation = "peak_per_orange_vs_footprint_length/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_per_methylation_for_site_{site}_sam_flag_{sam_flag}.tsv", # per: percentage 
        read_per_methylation_pkl = "peak_per_orange_vs_footprint_length/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_per_methylation_for_site_{site}_sam_flag_{sam_flag}.pkl", # per: percentage 
        
    shell:
        "python scripts/get_real_footprint_lengths_for_not_clustered_reads.py"
        " {input.peak_data} {input.methylation_matrix_real_footprint_len_pkl}"
        " {output.read_per_orange_and_length} {output.read_per_methylation}" 
        " {output.read_per_methylation_pkl}"
        " {input.footprint_capped_pkl}"
        " {params.neighbor_left_from_peak_center} {params.neighbor_right_from_peak_center}"
        " {params.strand} {wildcards.lflank} {wildcards.rflank}"

rule read_based_binding_assignment: 
    input:
        read_per_orange_and_length = "peak_per_orange_vs_footprint_length/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_po_vs_footplen_for_site_{site}_sam_flag_{sam_flag}.tsv", # po: percentage orange; footplen: footprint length 
        read_per_methylation_pkl = "peak_per_orange_vs_footprint_length/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_per_methylation_for_site_{site}_sam_flag_{sam_flag}.pkl", # per: percentage 
        peak_data = "peak_footprints/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_sam_flag_{sam_flag}.tsv", 
        occluded_pkl = "occluded_edges_on_methylation_vec/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_min_flen_10_wobble_gap_1_occluded.pkl", 
        footprint_capped_pkl = "footprint_min_size_selection/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_tf_footprint_min_length_10_wobble_gap_1.pkl", 
        
        
    output:
        binding_labeled_reads_cond1 = "binding_labeled_reads_for_peak/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_binding_label.cond1.tsv",
        binding_labeled_reads_cond2 = "binding_labeled_reads_for_peak/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_binding_label.cond2.tsv",
        binding_labeled_reads_cond3 = "binding_labeled_reads_for_peak/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_binding_label.cond3.tsv",
        binding_labeled_reads_cond4 = "binding_labeled_reads_for_peak/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_binding_label.cond4.tsv",
        binding_labeled_reads_cond5 = "binding_labeled_reads_for_peak/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_binding_label.cond5.tsv",
        binding_labeled_reads_cond6 = "binding_labeled_reads_for_peak/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_binding_label.cond6.tsv",
        binding_labeled_reads_cond7 = "binding_labeled_reads_for_peak/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_binding_label.cond7.tsv",
        binding_labeled_reads_cond8 = "binding_labeled_reads_for_peak/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_binding_label.cond8.tsv",
        
    shell:
        #"python scripts/assign_binding_for_each_read.py"
        "python scripts/assign_binding_for_each_read_multiple_conditions.py"
        " {input.read_per_orange_and_length} {input.peak_data}" 
        " {output.binding_labeled_reads_cond1}"
        " {output.binding_labeled_reads_cond2}"
        " {output.binding_labeled_reads_cond3}"
        " {output.binding_labeled_reads_cond4}"
        " {output.binding_labeled_reads_cond5}"
        " {output.binding_labeled_reads_cond6}" 
        " {input.read_per_methylation_pkl}" 
        " {output.binding_labeled_reads_cond7}" 
        " {output.binding_labeled_reads_cond8}"
        " {input.occluded_pkl} {input.footprint_capped_pkl}"
rule extend_read_after_binding_assignment: 
    input:
        binding_labeled_reads = "binding_labeled_reads_for_peak/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_binding_label.{condition}.tsv", 
        footprint_capped_pkl = "footprint_min_size_selection/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_tf_footprint_min_length_10_wobble_gap_1.pkl"
    
    params:
        strand = lambda wildcards: get_strand_info(wildcards)
    output: 
        extended_binding_labeled_reads = "extended_binding_labeled_reads_for_peak/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_extend_from_peak_left_{lextend}_right_{rextend}.{condition}.tsv",
        extended_binding_labeled_reads_verbose = "extended_binding_labeled_reads_for_peak/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_extend_from_peak_left_{lextend}_right_{rextend}.{condition}.verbse.tsv"
    shell:
        "python scripts/extend_footprint.py {input.binding_labeled_reads}"
        " {input.footprint_capped_pkl} {wildcards.lflank} {wildcards.rflank}"
        " {wildcards.lextend} {wildcards.rextend}"
        " {output.extended_binding_labeled_reads}"
        " {params.strand}"
        " {output.extended_binding_labeled_reads_verbose}" 

rule plot_mnase_and_extended_footprints: 
    input:
        extended_binding_labeled_reads = "extended_binding_labeled_reads_for_peak/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_extend_from_peak_left_{lextend}_right_{rextend}.{condition}.tsv",
        site_mnase_data = "site_specific_mnase_data/taken_from_{ref_flank}_cluster_{cl_id}_span_from_center_{span}_vline_left_{left}_right_{right}_site_{site}.tsv",
    params:
        title =  lambda wildcards: "MNase Cluster " + wildcards.cl_id + " " +  config["site_annotation"][wildcards.cl_id][wildcards.site], 
        ylabel = lambda wildcards: config ["sam_flag_annotation"][wildcards.sam_flag], 
        strand = lambda wildcards: get_strand_info(wildcards)
    output:
        as_circles = "plots/read_based_clustering/as_circles_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_extend_from_peak_left_{lextend}_right_{rextend}_mnase_ref_{ref_flank}_span_{span}_vline_left_l_{left}_r_{right}.{condition}.extend.pdf",
        average_in_one = "plots/read_based_clustering/average_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_extend_from_peak_left_{lextend}_right_{rextend}_mnase_ref_{ref_flank}_span_{span}_vline_left_l_{left}_r_{right}.{condition}.extend.pdf",
        averge_tsv = "plots/read_based_clustering/average_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_extend_from_peak_left_{lextend}_right_{rextend}_mnase_ref_{ref_flank}_span_{span}_vline_left_l_{left}_r_{right}.{condition}.extend.tsv",
    shell:
        "Rscript scripts/draw_methylation_status_v3_extend.R"
        " {input.extended_binding_labeled_reads} {output.as_circles}"
        " {wildcards.lextend} {wildcards.rextend} \"{params.title}\""
        " \"{params.ylabel}\" {output.average_in_one} {output.averge_tsv}"
        " \"{params.strand}\" {input.site_mnase_data}"
        " {wildcards.lflank} {wildcards.rflank}"
        " {wildcards.site} {wildcards.sam_flag}"

rule data_for_scatter_plot_by_binding: 
    input:
        binding_labeled_reads = "binding_labeled_reads_for_peak/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_binding_label.{condition}.tsv", 
        methylation_matrix_real_footprint_len_pkl =\
              "flank_footprint_matrix/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_methylation_matrix_flen.pkl",
    params:
    output:
        binding_labelled_scatter_tsv = "scatter_plot_by_binding/scatter_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}.{condition}.tsv"
    shell:
        "python scripts/scatter_plot_by_binding.py {input.binding_labeled_reads}"
        " {input.methylation_matrix_real_footprint_len_pkl}"
        " {output.binding_labelled_scatter_tsv}"

rule scatter_and_hist_plot_by_binding:
    input:
        binding_labelled_scatter_tsv = "scatter_plot_by_binding/scatter_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}.{condition}.tsv"
    output:
        scatter_plot_png = "plots/scatter_plot_by_binding/scatter_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}.{condition}.png", 
        hist_footprint_png = "plots/scatter_plot_by_binding/hist_footprint_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}.{condition}.png", 
        hist_orange_png = "plots/scatter_plot_by_binding/hist_orange_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}.{condition}.png", 
    shell:
        "Rscript scripts/scatter_plot_by_binding.R {input.binding_labelled_scatter_tsv}"
        " {output.scatter_plot_png} {output.hist_footprint_png} {output.hist_orange_png}"

rule combined_plots_read_based_binding:
    input:
        as_circles_99 = "plots/read_based_clustering/as_circles_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_99~147_extend_from_peak_left_{lextend}_right_{rextend}_mnase_ref_{ref_flank}_span_{span}_vline_left_l_{left}_r_{right}.{condition}.extend.pdf",
        as_circles_83 = "plots/read_based_clustering/as_circles_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_83~163_extend_from_peak_left_{lextend}_right_{rextend}_mnase_ref_{ref_flank}_span_{span}_vline_left_l_{left}_r_{right}.{condition}.extend.pdf",
        scatter_plot_png_83 = "plots/scatter_plot_by_binding/scatter_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_83~163.{condition}.png",
        scatter_plot_png_99 = "plots/scatter_plot_by_binding/scatter_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_99~147.{condition}.png",
        hist_footprint_png_83 = "plots/scatter_plot_by_binding/hist_footprint_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_83~163.{condition}.png", 
        hist_footprint_png_99 = "plots/scatter_plot_by_binding/hist_footprint_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_99~147.{condition}.png",
        hist_orange_png_83 = "plots/scatter_plot_by_binding/hist_orange_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_83~163.{condition}.png", 
        hist_orange_png_99 = "plots/scatter_plot_by_binding/hist_orange_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_99~147.{condition}.png",
        
        
    output:
         combined_png =  "plots/combined_plots_read_based_clustering/as_circles_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_extend_from_peak_left_{lextend}_right_{rextend}_mnase_ref_{ref_flank}_span_{span}_vline_l_{left}_r_{right}.combined.{condition}.png"
    shell:
        "sh scripts/combined_read_based_plots.sh"
        " {input.as_circles_99} {input.as_circles_83}"
        " {input.scatter_plot_png_99} {input.scatter_plot_png_83}"
        " {input.hist_footprint_png_99} {input.hist_footprint_png_83}"
        " {input.hist_orange_png_99} {input.hist_orange_png_83}" 
        " {output.combined_png}" 





################## Generating V-plot Cluster 3 matrix ####### 
def get_files_for_occupancy_matrix(wildcards):
    flist = [] 
    peak_fp = open(config[wildcards.setting])
    for line in peak_fp:
        line_items = line.strip().split()
        peak_name = line_items[0]
        cl_id = line_items[1] 
        fname = "binding_labeled_reads_for_peak/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_binding_label.{condition}.tsv".format(
           sample = wildcards.sample, 
           bed = wildcards.bed,
           cl_id = cl_id,
           lflank = wildcards.lflank,
           rflank = wildcards.rflank, 
           site = peak_name, 
           sam_flag = wildcards.sam_flag,
           condition = wildcards.condition) 
        flist.append(fname)
    return flist
rule get_file_list_for_occupancy_matrix: 
    input:
        all_files = lambda wildcards: get_files_for_occupancy_matrix (wildcards)
    output:
        flist_occupancy_matrix = "flist_occupancy_matrix_for_vplot_clusters/{sample}_to_{bed}_occup_matrix_for_{setting}_with_lf_{lflank}_rf_{rflank}_sam_flag_{sam_flag}_condition_{condition}.flist.tsv"
    run:
        fp = open (output.flist_occupancy_matrix, "w")
        for f in input.all_files:
            fp.write(f + "\n")
        fp.close()


rule build_occupancy_matrix: 
    input:
        flist_occupancy_matrix = "flist_occupancy_matrix_for_vplot_clusters/{sample}_to_{bed}_occup_matrix_for_{setting}_with_lf_{lflank}_rf_{rflank}_sam_flag_{sam_flag}_condition_{condition}.flist.tsv"
    
    output:
        occupancy_matrix = "occupancy_matrix_for_vplot_clusters/{sample}_to_{bed}_occup_matrix_for_{setting}_with_lf_{lflank}_rf_{rflank}_sam_flag_{sam_flag}_condition_{condition}.tsv"
        
    shell:
        "python scripts/preapre_occupancy_matrix_vplot_cl3.py"
        " {input.flist_occupancy_matrix} {output.occupancy_matrix}"


rule build_merged_occupancy_matrix:
    input:
        flist_occupancy_matrix_83 = "flist_occupancy_matrix_for_vplot_clusters/{sample}_to_{bed}_occup_matrix_for_{setting}_with_lf_{lflank}_rf_{rflank}_sam_flag_83~163_condition_{condition}.flist.tsv", 
        flist_occupancy_matrix_99 = "flist_occupancy_matrix_for_vplot_clusters/{sample}_to_{bed}_occup_matrix_for_{setting}_with_lf_{lflank}_rf_{rflank}_sam_flag_99~147_condition_{condition}.flist.tsv"
    output:
        merged_occupancy_matrix = "occupancy_matrix_for_vplot_clusters/{sample}_to_{bed}_occup_matrix_for_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.merged.tsv"
    shell:
        "python scripts/preapre_merged_occupancy_matrix_vplot_cl3.py"
        " {input.flist_occupancy_matrix_83} {input.flist_occupancy_matrix_99}" 
        " {output.merged_occupancy_matrix}"

rule kmeans_on_meged_occupancy_matrix:
    input:
        merged_occupancy_matrix = "occupancy_matrix_for_vplot_clusters/{sample}_to_{bed}_occup_matrix_for_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.merged.tsv"
    output:
        kmeans_tsv = "kmeans_on_occupany/{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.tsv",
        kmeans_center = "kmeans_on_occupany/{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.center.tsv",
        kmeans_cluster_cnt = "kmeans_on_occupany/{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.cnt.tsv",
        kmeans_rdata = "kmeans_on_occupany/{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.Rdata", 
        kmeans_original_complete = "kmeans_on_occupany/{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.ori.merged.tsv", 
        kmeans_original_center = "kmeans_on_occupany/{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.ori.center.tsv",
        kmeans_original_cluster_cnt = "kmeans_on_occupany/{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.ori.cnt.tsv",

    shell:
        "Rscript scripts/kmeans_on_occupancy_matrix.R"
        " {input.merged_occupancy_matrix} {wildcards.nclust} {output.kmeans_tsv}"
        " {output.kmeans_rdata} {output.kmeans_center}"
        " {output.kmeans_cluster_cnt}"
        " {output.kmeans_original_complete} {output.kmeans_original_center}"
        " {output.kmeans_original_cluster_cnt}" 
    
rule center_occupancy_kmeans:
    input:
        kmeans_center = "kmeans_on_occupany/{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.center.tsv",
        kmeans_cluster_cnt = "kmeans_on_occupany/{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.cnt.tsv",
        gnuplt_base_file = "utils/gnuplot_base_files/center_kmeans.gplt"    
    params:
        plt_title = lambda wildcards: config["plot_titles"][wildcards.setting] + " (" + wildcards.condition + ")" 
    output:
        center_gplt = "plots/kmeans_on_occupancy/{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.center.gplt",
        center_eps = "plots/kmeans_on_occupancy/{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.center.eps",
        center_pdf = "plots/kmeans_on_occupancy/{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.center.pdf",
    shell:
        "sh scripts/kmeans_center_plot.sh {input.kmeans_center}"
        " {input.gnuplt_base_file} {output.center_gplt} {output.center_eps}"
        " {output.center_pdf} \"{params.plt_title}\" {input.kmeans_cluster_cnt}"


rule plot_whole_occupancy_matrix: 
    input:
        kmeans_tsv = "kmeans_on_occupany/{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.tsv",
        kmeans_cluster_cnt = "kmeans_on_occupany/{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.cnt.tsv",
        gnuplt_base_file = "utils/gnuplot_base_files/occupancy_heatmap.gplt"
    params:
        plt_title = lambda wildcards: config["plot_titles"][wildcards.setting] + " (" + wildcards.condition + ")" 
    output:
        occupancy_heatmap_gplt = "plots/kmeans_on_occupancy/{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.complete.gplt" , 
        occupancy_heatmap_eps = "plots/kmeans_on_occupancy/{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.complete.eps",
        occupancy_heatmap_pdf = "plots/kmeans_on_occupancy/{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.complete.pdf",
        
    shell:
        "sh scripts/kmeans_complete_data_heatmap.sh {input.kmeans_tsv}"
        " {input.gnuplt_base_file} {output.occupancy_heatmap_gplt} {output.occupancy_heatmap_eps}"
        " {output.occupancy_heatmap_pdf} \"{params.plt_title}\" {input.kmeans_cluster_cnt}"

rule plot_whole_occupancy_matrix_no_sorting:
    input:
        kmeans_tsv = "kmeans_on_occupany/{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.ori.merged.tsv",
        kmeans_cluster_cnt = "kmeans_on_occupany/{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.ori.cnt.tsv",
        gnuplt_base_file = "utils/gnuplot_base_files/occupancy_heatmap.gplt"
    params:
        plt_title = lambda wildcards: config["plot_titles"][wildcards.setting] + " (" + wildcards.condition + ")" 
    output:
        occupancy_heatmap_gplt = "plots/kmeans_on_occupancy/{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.ori.complete.gplt" , 
        occupancy_heatmap_eps = "plots/kmeans_on_occupancy/{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.ori.complete.eps",
        occupancy_heatmap_pdf = "plots/kmeans_on_occupancy/{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.ori.complete.pdf",
        
    shell:
        "sh scripts/kmeans_complete_data_heatmap.sh {input.kmeans_tsv}"
        " {input.gnuplt_base_file} {output.occupancy_heatmap_gplt} {output.occupancy_heatmap_eps}"
        " {output.occupancy_heatmap_pdf} \"{params.plt_title}\" {input.kmeans_cluster_cnt}"

rule center_occupancy_kmeans_no_sorting:
    input:
        kmeans_center = "kmeans_on_occupany/{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.ori.center.tsv",
        kmeans_cluster_cnt = "kmeans_on_occupany/{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.ori.cnt.tsv",
        gnuplt_base_file = "utils/gnuplot_base_files/center_kmeans.gplt"    
    params:
        plt_title = lambda wildcards: config["plot_titles"][wildcards.setting] + " (" + wildcards.condition + ")" 
    output:
        center_gplt = "plots/kmeans_on_occupancy/{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.ori.center.gplt",
        center_eps = "plots/kmeans_on_occupancy/{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.ori.center.eps",
        center_pdf = "plots/kmeans_on_occupancy/{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.ori.center.pdf",
    shell:
        "sh scripts/kmeans_center_plot.sh {input.kmeans_center}"
        " {input.gnuplt_base_file} {output.center_gplt} {output.center_eps}"
        " {output.center_pdf} \"{params.plt_title}\" {input.kmeans_cluster_cnt}"


    
def get_files_for_elbow_for_occupancy(wildcards):
    flist = []
    for cl in config["kmeans_settings"][wildcards.setting]:
        flist.append("kmeans_on_occupany/{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.Rdata".format(
                      sample = wildcards.sample,
                      bed = wildcards.bed,
                      setting = wildcards.setting, 
                      lflank = wildcards.lflank,
                      rflank = wildcards.rflank,
                      condition = wildcards.condition, 
                      nclust = cl))
    return flist 
rule elbow_occupancy_kmeans: 
    input:
        rdata_files = lambda wildcards: get_files_for_elbow_for_occupancy(wildcards) 
    params:
        cl_ids = lambda wildcards: "@".join(map(str,
                     config["kmeans_settings"][wildcards.setting]))
    output:
        elbow_tsv = "plots/kmeans_on_occupancy/{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.elbow.tsv",
        elbow_png = "plots/kmeans_on_occupancy/{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.elbow.png",
        
    shell:
        "Rscript scripts/elbow_chip_kmeans.R \"{input.rdata_files}\" {params.cl_ids} {output.elbow_tsv}"
        " {output.elbow_png} {wildcards.setting}"


rule prepare_bed_for_mnase_plot:
    input:
        kmeans_tsv = "kmeans_on_occupany/{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.tsv",
        input_bed = "input_bed/reordered_20clust_147_50PNE_open_v1.bed",
        genome_size_file = "metadata/dm3.chrom.sizes"
    params:
    output:
        bed_for_mnase = "mnase_based_on_occupancy/{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.slop_{slop}.bed"
        
    shell:
        "sh scripts/selected_for_peaks.sh {input.kmeans_tsv} {input.input_bed} {output.bed_for_mnase} {wildcards.slop} {input.genome_size_file}"


rule map_manse_bigwig_to_bed:
    input:
        bed_for_mnase = "mnase_based_on_occupancy/{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.slop_{slop}.bed", 
        bw_file = "bigwigs/{bw}.bigwig",
    params:
        strand_column = 6
    output:
        raw_csv_gz = "mnase_based_on_occupancy/{bw}_mapped_on_{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.slop_{slop}.csv.gz", 
        eom_csv_gz = "mnase_based_on_occupancy/{bw}_mapped_on_{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.slop_{slop}.csv.gz_e.csv.gz",
        exat_stats_tsv = "mnase_based_on_occupancy/{bw}_mapped_on_{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.slop_{slop}.csv.gz_exact_stats.tsv"
    shell:
        "python $NGS_SCRIPTS_DIR/map_bw_to_bed_strand_aware.py"
        " {input.bw_file} {input.bed_for_mnase} {output.raw_csv_gz} {params.strand_column}" 
    

rule mean_mnase_plots_for_clusters:
    input:
        eom_csv_gz = "mnase_based_on_occupancy/{bw}_mapped_on_{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.slop_{slop}.csv.gz_e.csv.gz",
    params:
    output:
        colmeans_tsv = "mnase_based_on_occupancy/{bw}_mapped_on_{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.slop_{slop}.colmeans.tsv",  
    shell:
        "sh scripts/calculate_colmeans_of_mnase_eom.sh {input.eom_csv_gz} {wildcards.nclust} {wildcards.slop} {output.colmeans_tsv} " 

rule plot_mean_mnase_for_occupancy:
    input:
        colmeans_tsv = "mnase_based_on_occupancy/{bw}_mapped_on_{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.slop_{slop}.colmeans.tsv", 
        base_gnuplt_file = "utils/gnuplot_base_files/mnase_colmeans.gplt", 
        
    params:
        cluster_id_column = 3 
    output:
        mnase_colmean_plot_pdf = "plots/mnase_map_on_occupancy/{bw}_mapped_on_{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.slop_{slop}.colmeans.pdf", 
        mnase_colmean_plot_png = "plots/mnase_map_on_occupancy/{bw}_mapped_on_{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.slop_{slop}.colmeans.png", 
        mnase_colmean_gplt = "plots/mnase_map_on_occupancy/{bw}_mapped_on_{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.slop_{slop}.colmeans.gplt", 

    shell:
        #"Rscript scripts/colmeans_mnase_occupancy_plot.R {input.colmeans_tsv}" # Rcode  
        #" {output.mnase_colmean_plot_pdf} {output.mnase_colmean_plot_png}"
        #" {wildcards.bw}" 
        "sh scripts/colmeans_mnase_occupancy_plot.sh {input.colmeans_tsv}"
        " {params.cluster_id_column} {output.mnase_colmean_plot_pdf} {output.mnase_colmean_plot_png}  {input.base_gnuplt_file} {output.mnase_colmean_gplt} "
        " {wildcards.bw} {wildcards.slop}" 


rule short_vs_long_mnase_for_occupancy:
    input:
        short_mnase_colmean = "mnase_based_on_occupancy/{short_bw}_mapped_on_{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.slop_{slop}.colmeans.tsv",
        long_mnase_colmean = "mnase_based_on_occupancy/{long_bw}_mapped_on_{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.slop_{slop}.colmeans.tsv",
        kmeans_cluster_cnt = "kmeans_on_occupany/{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.cnt.tsv", 
        base_gnuplt_file = "utils/gnuplot_base_files/mnase_colmeans.gplt", 

    params:
        plot_title = lambda wildcards: "Cluster " + wildcards.kclust_id + ": " + wildcards.short_bw + " vs. " + wildcards.long_bw 
    output:
        short_vs_long_colmeans_pdf = "plots/mnase_map_on_occupancy/{short_bw}_vs_{long_bw}_mapped_on_{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.slop_{slop}.kclust_{kclust_id}_colmeans.pdf",
        short_vs_long_colmeans_png = "plots/mnase_map_on_occupancy/{short_bw}_vs_{long_bw}_mapped_on_{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.slop_{slop}.kclust_{kclust_id}_colmeans.png",
        short_vs_long_colmeans_gplt = "plots/mnase_map_on_occupancy/{short_bw}_vs_{long_bw}_mapped_on_{sample}_to_{bed}_ocm_{setting}_with_lf_{lflank}_rf_{rflank}_condition_{condition}.kmeans_nclust_{nclust}.merged.slop_{slop}.kclust_{kclust_id}_colmeans.gplt",
    shell:
        # "Rscript scripts/colmeans_mnase_combined_occupancy_plot.R"
        "sh scripts/colmeans_mnase_combined_occupancy_plot.sh"
        " {input.short_mnase_colmean} {input.long_mnase_colmean}"
        " {input.kmeans_cluster_cnt} {input.base_gnuplt_file}"
        " {output.short_vs_long_colmeans_gplt} {output.short_vs_long_colmeans_pdf}"
        " {output.short_vs_long_colmeans_png} {wildcards.kclust_id} "
        " \"{params.plot_title}\"" 

