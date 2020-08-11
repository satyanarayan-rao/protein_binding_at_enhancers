######### Perform Clustering using footprint vector at bp resolution ####### 
rule get_footprint_length_vec_per_bp_for_peak:
    input:
        footprint_length_at_bp_resolution_tsv = \
              "flank_footprint_matrix/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_methylation_matrix.footprint_len.bp.res.tsv", 
    params:
        string_for_site = lambda wildcards: config["site_annotation"][wildcards.cl_id][wildcards.site] + "\^" + wildcards.cl_id
    output:
        peak_footprint_length_at_bp_resolution_tsv = "peak_specific_footprint_length_vector/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_footprint_length_matrix_site_{site}_sam_flag_{sam_flag}.tsv"
    shell:
        "grep -w \"{params.string_for_site}\" {input.footprint_length_at_bp_resolution_tsv} | grep -w \"{wildcards.sam_flag}\" > {output.peak_footprint_length_at_bp_resolution_tsv}"

rule cluster_footprint_length_vec:
    input:
        peak_footprint_length_at_bp_resolution_tsv = "peak_specific_footprint_length_vector/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_footprint_length_matrix_site_{site}_sam_flag_{sam_flag}.tsv",
    params:
        consider_from_center = 10
    output:
        clustered_vec = "kmeans_peak_specific_footprint_length_vector/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_footprint_length_matrix_site_{site}_sam_flag_{sam_flag}_nclust_{nclust}.tsv", 
        read_to_cl = "kmeans_peak_specific_footprint_length_vector/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_footprint_length_matrix_site_{site}_sam_flag_{sam_flag}_nclust_{nclust}.read.to.label_map.tsv", 
        kmeans_rdata = "kmeans_peak_specific_footprint_length_vector/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_footprint_length_matrix_site_{site}_sam_flag_{sam_flag}_nclust_{nclust}.read.to.label_map.Rdata",
    shell:
        "Rscript scripts/kmeans_footprint_length_vec.R"
        " {input.peak_footprint_length_at_bp_resolution_tsv}"
        " {output.clustered_vec} {output.read_to_cl}"
        " {wildcards.lflank} {wildcards.rflank} {params.consider_from_center}" 
        " {wildcards.nclust} {output.kmeans_rdata}"

def get_kmeans_rdata_files (wildcards):
    flist = []
    for cl in config["kmeans_settings"][wildcards.setting]:
        flist.append("kmeans_peak_specific_footprint_length_vector/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_footprint_length_matrix_site_{site}_sam_flag_{sam_flag}_nclust_{nclust}.read.to.label_map.Rdata".format(
                      sample = wildcards.sample,
                      bed = wildcards.bed,
                      cl_id = wildcards.cl_id,
                      lflank = wildcards.lflank,
                      rflank = wildcards.rflank,
                      site = wildcards.site, 
                      sam_flag = wildcards.sam_flag,
                      nclust = cl))
    return flist

rule elbow_analyis_chip_scores:
    input:
        rdata_files = lambda wildcards: get_kmeans_rdata_files (wildcards)
    params:
        cl_ids = lambda wildcards: "@".join(map(str,
                     config["kmeans_settings"][wildcards.setting]))
    output:
        elbow_tsv = "plots/kmeans_footprint/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_footprint_length_matrix_site_{site}_sam_flag_{sam_flag}_setting_{setting}.elbow.tsv",
        elbow_png = "plots/kmeans_footprint/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_footprint_length_matrix_site_{site}_sam_flag_{sam_flag}_setting_{setting}.elbow.png"
    shell:
        "Rscript scripts/elbow_chip_kmeans.R \"{input.rdata_files}\" {params.cl_ids} {output.elbow_tsv}"
        " {output.elbow_png} {wildcards.setting}"

rule annotate_footprint_vec_with_kmeans: 
    input:
        peak_data = "peak_footprints/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_sam_flag_{sam_flag}.tsv",  
        read_to_cl = "kmeans_peak_specific_footprint_length_vector/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_footprint_length_matrix_site_{site}_sam_flag_{sam_flag}_nclust_{nclust}.read.to.label_map.tsv" 
    params:
    output:
        annotated_vec = "cl_annotated_peak_footprints/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_sam_flag_{sam_flag}_nclust_{nclust}.with_label.tsv"
    shell:
        "python scripts/assign_kmeans_label.py {input.read_to_cl} {input.peak_data}"
        " {output.annotated_vec}"

rule extend_kmeans_annotated_vector:
    input:
        annotated_vec = "cl_annotated_peak_footprints/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_sam_flag_{sam_flag}_nclust_{nclust}.with_label.tsv",
        extended_matrix_pkl = "footprint_min_size_selection/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_tf_footprint_min_length_10_wobble_gap_1.pkl"
    params:
        strand = lambda wildcards: get_strand_info(wildcards)
    output:
        extended_footprint_tsv = "extended_kmeans_footprints/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_sam_flag_{sam_flag}_nclust_{nclust}_extend_from_peak_left_{lextend}_right_{rextend}.tsv" 
    shell:
        "python scripts/extend_footprint.py {input.annotated_vec}"
        " {input.extended_matrix_pkl} {wildcards.lflank} {wildcards.rflank}"
        " {wildcards.lextend} {wildcards.rextend} {output.extended_footprint_tsv}"
        " {params.strand}"

rule data_for_scatter_kmeans:
    input:
        annotated_vec = "cl_annotated_peak_footprints/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_sam_flag_{sam_flag}_nclust_{nclust}.with_label.tsv",
        methylation_matrix_real_footprint_len_pkl =\
              "flank_footprint_matrix/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_methylation_matrix_flen.pkl",
        
    params:
    output:
        binding_labelled_scatter_tsv = "scatter_plot_by_kmeans/scatter_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_nclust_{nclust}.tsv"
    shell:
        "python scripts/scatter_plot_by_binding.py {input.annotated_vec}"
        " {input.methylation_matrix_real_footprint_len_pkl}"
        " {output.binding_labelled_scatter_tsv}" 

rule scatter_and_hist_plot_kmeans:
    input:
        binding_labelled_scatter_tsv = "scatter_plot_by_kmeans/scatter_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_nclust_{nclust}.tsv"
    output:
        scatter_plot_png = "plots/scatter_plot_by_kmeans/scatter_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_nclust_{nclust}.png", 
        hist_footprint_png = "plots/scatter_plot_by_kmeans/hist_footprint_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_nclust_{nclust}.png", 
        hist_orange_png = "plots/scatter_plot_by_kmeans/hist_orange_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_nclust_{nclust}.png", 
    shell:
        "Rscript scripts/scatter_plot_by_kmeans.R {input.binding_labelled_scatter_tsv}"
        " {output.scatter_plot_png} {output.hist_footprint_png} {output.hist_orange_png}"

rule plot_mnase_and_extended_footprints_kmeans: 
    input:
        extended_footprint_tsv = "extended_kmeans_footprints/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_sam_flag_{sam_flag}_nclust_{nclust}_extend_from_peak_left_{lextend}_right_{rextend}.tsv",
        site_mnase_data = "site_specific_mnase_data/taken_from_{ref_flank}_cluster_{cl_id}_span_from_center_{span}_vline_left_{left}_right_{right}_site_{site}.tsv",
    params:
        title =  lambda wildcards: "MNase Cluster " + wildcards.cl_id + " " +  config["site_annotation"][wildcards.cl_id][wildcards.site], 
        ylabel = lambda wildcards: config ["sam_flag_annotation"][wildcards.sam_flag], 
        strand = lambda wildcards: get_strand_info(wildcards)
    output:
        as_circles = "plots/kmeans_based_clustering/as_circles_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_nclust_{nclust}_extend_from_peak_left_{lextend}_right_{rextend}_mnase_ref_{ref_flank}_span_{span}_vline_left_l_{left}_r_{right}.extend.pdf",
        average_in_one = "plots/kmeans_based_clustering/average_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_nclust_{nclust}_extend_from_peak_left_{lextend}_right_{rextend}_mnase_ref_{ref_flank}_span_{span}_vline_left_l_{left}_r_{right}.extend.pdf",
        averge_tsv = "plots/kmeans_based_clustering/average_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_nclust_{nclust}_extend_from_peak_left_{lextend}_right_{rextend}_mnase_ref_{ref_flank}_span_{span}_vline_left_l_{left}_r_{right}.extend.tsv",
    shell:
        "Rscript scripts/draw_methylation_status_v3_extend.R"
        " {input.extended_footprint_tsv} {output.as_circles}"
        " {wildcards.lextend} {wildcards.rextend} \"{params.title}\""
        " \"{params.ylabel}\" {output.average_in_one} {output.averge_tsv}"
        " \"{params.strand}\" {input.site_mnase_data}"
        " {wildcards.lflank} {wildcards.rflank}"
        " {wildcards.site} {wildcards.sam_flag}"

rule combined_plots_read_based_kmeans:
    input:
        as_circles_99 = "plots/kmeans_based_clustering/as_circles_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_99~147_nclust_{nclust}_extend_from_peak_left_{lextend}_right_{rextend}_mnase_ref_{ref_flank}_span_{span}_vline_left_l_{left}_r_{right}.extend.pdf",
        as_circles_83 = "plots/kmeans_based_clustering/as_circles_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_83~163_nclust_{nclust}_extend_from_peak_left_{lextend}_right_{rextend}_mnase_ref_{ref_flank}_span_{span}_vline_left_l_{left}_r_{right}.extend.pdf",
        scatter_plot_png_83 = "plots/scatter_plot_by_kmeans/scatter_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_83~163_nclust_{nclust}.png",
        scatter_plot_png_99 = "plots/scatter_plot_by_kmeans/scatter_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_99~147_nclust_{nclust}.png",
        hist_footprint_png_83 = "plots/scatter_plot_by_kmeans/hist_footprint_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_83~163_nclust_{nclust}.png", 
        hist_footprint_png_99 = "plots/scatter_plot_by_kmeans/hist_footprint_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_99~147_nclust_{nclust}.png",
        hist_orange_png_83 = "plots/scatter_plot_by_kmeans/hist_orange_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_83~163_nclust_{nclust}.png", 
        hist_orange_png_99 = "plots/scatter_plot_by_kmeans/hist_orange_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_99~147_nclust_{nclust}.png",
        
    output:
         combined_png =  "plots/combined_plots_kmeans/as_circles_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_nclust_{nclust}_extend_from_peak_left_{lextend}_right_{rextend}_mnase_ref_{ref_flank}_span_{span}_vline_l_{left}_r_{right}.combined.png"
    shell:
        "sh scripts/combined_read_based_plots.sh"
        " {input.as_circles_99} {input.as_circles_83}"
        " {input.scatter_plot_png_99} {input.scatter_plot_png_83}"
        " {input.hist_footprint_png_99} {input.hist_footprint_png_83}"
        " {input.hist_orange_png_99} {input.hist_orange_png_83}" 
        " {output.combined_png}" 
