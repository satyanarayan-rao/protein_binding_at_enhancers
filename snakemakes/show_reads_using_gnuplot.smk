#rule put_e_for_the_edges_on_extended_vectors:
#     input:
#         extended_binding_labeled_reads_verbose = "extended_binding_labeled_reads_for_peak/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_extend_from_peak_left_{lextend}_right_{rextend}.{condition}.verbse.tsv",
#         occluded_pkl = "occluded_edges_on_methylation_vec/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_min_flen_10_wobble_gap_1_occluded.pkl", 
#     output:
#         extended_binding_labeled_reads_verbose_edges_with_e = "extended_binding_labeled_reads_for_peak/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_extend_from_peak_left_{lextend}_right_{rextend}.{condition}.edges_verbse.tsv", 
#     shell:
#         "python scripts/add_e_at_edges.py"
#         " {input.extended_binding_labeled_reads_verbose}"
#         " {input.occluded_pkl}"
#         " {output.extended_binding_labeled_reads_verbose_edges_with_e}"

rule order_extended_footprint_vec_single_binding_to_matrix:
    input:
        extended_binding_labeled_reads_verbose = "extended_binding_labeled_reads_for_peak/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_extend_from_peak_left_{lextend}_right_{rextend}.{condition}.verbse.tsv"
        #extended_binding_labeled_reads_verbose_edges_with_e = "extended_binding_labeled_reads_for_peak/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_extend_from_peak_left_{lextend}_right_{rextend}.{condition}.edges_verbse.tsv", 
    params:
        
    output:
        ordered_footprint_vec = "extended_binding_labeled_reads_for_peak/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_extend_from_peak_left_{lextend}_right_{rextend}.{condition}.footprint_ordered.tsv",
        ordered_methylation_vec = "extended_binding_labeled_reads_for_peak/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_extend_from_peak_left_{lextend}_right_{rextend}.{condition}.mvec_ordered.tsv"
    shell:
        "sh scripts/order_methylation_and_footprint_vec_single_binding.sh"
        #" {input.extended_binding_labeled_reads_verbose_edges_with_e}" 
        " {input.extended_binding_labeled_reads_verbose}" 
        " {output.ordered_footprint_vec} {output.ordered_methylation_vec}"

rule convert_dot_to_numeric_matrix_footprint:
    input:
        ordered_footprint_vec = "extended_binding_labeled_reads_for_peak/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_extend_from_peak_left_{lextend}_right_{rextend}.{condition}.footprint_ordered.tsv",
    output:
        ordered_footprint_matrix = "extended_binding_labeled_reads_for_peak/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_extend_from_peak_left_{lextend}_right_{rextend}.{condition}.footprint_num_matrix.tsv",

    shell:
        "python scripts/footprint_dot_to_digit_vec.py {input.ordered_footprint_vec}"
        " {output.ordered_footprint_matrix}" 

rule convert_dot_to_numeric_matrix_mvec:
    input:
        ordered_methylation_vec = "extended_binding_labeled_reads_for_peak/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_extend_from_peak_left_{lextend}_right_{rextend}.{condition}.mvec_ordered.tsv"
        
    output:
        ordered_methylation_matrix = "extended_binding_labeled_reads_for_peak/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_extend_from_peak_left_{lextend}_right_{rextend}.{condition}.methylation_num_matrix.tsv",
    shell:
        "python scripts/methylation_dot_to_digit_vec.py {input.ordered_methylation_vec}"
        " {output.ordered_methylation_matrix}"

rule plot_methylation_and_footprint_by_sam_flag:
    input:
        ordered_footprint_matrix = "extended_binding_labeled_reads_for_peak/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_extend_from_peak_left_{lextend}_right_{rextend}.{condition}.footprint_num_matrix.tsv",
        site_mnase_data = "site_specific_mnase_data/taken_from_{ref_flank}_cluster_{cl_id}_span_from_center_{span}_vline_left_{left}_right_{right}_site_{site}.tsv", 
        ordered_methylation_matrix = "extended_binding_labeled_reads_for_peak/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_extend_from_peak_left_{lextend}_right_{rextend}.{condition}.methylation_num_matrix.tsv", 

        gnuplt_mnase_params = "utils/gnuplot_base_files/manse_params.gplt", 
        gnuplt_footprint_params = "utils/gnuplot_base_files/footprint_params.gplt", 
        gnuplt_methylation_params = "utils/gnuplot_base_files/methylation_params.gplt", 
        
    output:
        ordered_footprint_pdf = "plots/extended_gnuplot_based_figures/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_extend_from_peak_left_{lextend}_right_{rextend}_mnase_ref_{ref_flank}_span_{span}_vline_left_l_{left}_r_{right}.{condition}.footprint_num_matrix.pdf", 
        ordered_methylation_pdf = "plots/extended_gnuplot_based_figures/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_extend_from_peak_left_{lextend}_right_{rextend}_mnase_ref_{ref_flank}_span_{span}_vline_left_l_{left}_r_{right}.{condition}.methylation_num_matrix.pdf", 
        ordered_footprint_gplt = "plots/extended_gnuplot_based_figures/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_extend_from_peak_left_{lextend}_right_{rextend}_mnase_ref_{ref_flank}_span_{span}_vline_left_l_{left}_r_{right}.{condition}.footprint_num_matrix.gplt", 
        ordered_methylation_gplt = "plots/extended_gnuplot_based_figures/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_{sam_flag}_extend_from_peak_left_{lextend}_right_{rextend}_mnase_ref_{ref_flank}_span_{span}_vline_left_l_{left}_r_{right}.{condition}.methylation_num_matrix.gplt", 
    shell: 
        "sh scripts/plot_footprint_and_methylation.sh"
        " {input.ordered_footprint_matrix} {input.ordered_methylation_matrix}"
        " {input.site_mnase_data} {input.gnuplt_mnase_params}" 
        " {input.gnuplt_footprint_params} {output.ordered_footprint_pdf}"
        " {output.ordered_methylation_pdf}"
        " {output.ordered_footprint_gplt} {output.ordered_methylation_gplt}"
        " {wildcards.lextend} {wildcards.rextend}"
        " {wildcards.lflank} {wildcards.rflank}"
        " {input.gnuplt_methylation_params}"


rule merge_sam_flags_footprint_and_methylated_vec:
    input:
        ordered_99_147_footprint = "extended_binding_labeled_reads_for_peak/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_99~147_extend_from_peak_left_{lextend}_right_{rextend}.{condition}.footprint_ordered.tsv",
        ordered_83_163_footprint = "extended_binding_labeled_reads_for_peak/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_83~163_extend_from_peak_left_{lextend}_right_{rextend}.{condition}.footprint_ordered.tsv", 
        ordered_99_147_methylation = "extended_binding_labeled_reads_for_peak/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_99~147_extend_from_peak_left_{lextend}_right_{rextend}.{condition}.mvec_ordered.tsv",
        ordered_83_163_methylation = "extended_binding_labeled_reads_for_peak/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_sam_flag_83~163_extend_from_peak_left_{lextend}_right_{rextend}.{condition}.mvec_ordered.tsv", 
    params:
        
    output:
        merged_footprints = "extended_binding_labeled_reads_for_peak/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_extend_from_peak_left_{lextend}_right_{rextend}.{condition}.footprint_merged_ordered.tsv",
        merged_methylation = "extended_binding_labeled_reads_for_peak/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_extend_from_peak_left_{lextend}_right_{rextend}.{condition}.mvec_merged_ordered.tsv" , 
         
    shell:
        "sh scripts/merge_and_order_footprints_and_methylation.sh"
        " {input.ordered_99_147_footprint} {input.ordered_83_163_footprint}"
        " {input.ordered_99_147_methylation} {input.ordered_83_163_methylation}"
        " {output.merged_footprints} {output.merged_methylation}" 

rule plot_merged_reads_footprint_and_methylation_together: 
    input:
        ordered_footprint_matrix = "extended_binding_labeled_reads_for_peak/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_extend_from_peak_left_{lextend}_right_{rextend}.{condition}.footprint_merged_ordered.tsv",
        site_mnase_data = "site_specific_mnase_data/taken_from_{ref_flank}_cluster_{cl_id}_span_from_center_{span}_vline_left_{left}_right_{right}_site_{site}.tsv", 
        ordered_methylation_matrix = "extended_binding_labeled_reads_for_peak/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_extend_from_peak_left_{lextend}_right_{rextend}.{condition}.mvec_merged_ordered.tsv", 

        gnuplt_mnase_params = "utils/gnuplot_base_files/manse_params.gplt", 
        gnuplt_footprint_params = "utils/gnuplot_base_files/footprint_params.gplt", 
        gnuplt_methylation_params = "utils/gnuplot_base_files/methylation_params.gplt", 
 
    params:
    output:
        ordered_footprint_pdf = "plots/extended_gnuplot_based_figures/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_extend_from_peak_left_{lextend}_right_{rextend}_mnase_ref_{ref_flank}_span_{span}_vline_left_l_{left}_r_{right}.{condition}.footprint_merged_num_matrix.pdf", 
        ordered_methylation_pdf = "plots/extended_gnuplot_based_figures/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_extend_from_peak_left_{lextend}_right_{rextend}_mnase_ref_{ref_flank}_span_{span}_vline_left_l_{left}_r_{right}.{condition}.methylation_merged_num_matrix.pdf", 
        ordered_footprint_gplt = "plots/extended_gnuplot_based_figures/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_extend_from_peak_left_{lextend}_right_{rextend}_mnase_ref_{ref_flank}_span_{span}_vline_left_l_{left}_r_{right}.{condition}.footprint_merged_num_matrix.gplt", 
        ordered_methylation_gplt = "plots/extended_gnuplot_based_figures/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_extend_from_peak_left_{lextend}_right_{rextend}_mnase_ref_{ref_flank}_span_{span}_vline_left_l_{left}_r_{right}.{condition}.methylation_merged_num_matrix.gplt",
        ordered_footprint_mat = "plots/extended_gnuplot_based_figures/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_extend_from_peak_left_{lextend}_right_{rextend}_mnase_ref_{ref_flank}_span_{span}_vline_left_l_{left}_r_{right}.{condition}.footprint_merged_num_matrix.mat.tsv", 
        ordered_methylation_mat = "plots/extended_gnuplot_based_figures/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_extend_from_peak_left_{lextend}_right_{rextend}_mnase_ref_{ref_flank}_span_{span}_vline_left_l_{left}_r_{right}.{condition}.methylation_merged_num_matrix.mat.tsv", 
        #ordered_footprint_eps = "plots/extended_gnuplot_based_figures/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_extend_from_peak_left_{lextend}_right_{rextend}_mnase_ref_{ref_flank}_span_{span}_vline_left_l_{left}_r_{right}.{condition}.footprint_merged_num_matrix.eps", 
        #ordered_methylation_eps = "plots/extended_gnuplot_based_figures/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_extend_from_peak_left_{lextend}_right_{rextend}_mnase_ref_{ref_flank}_span_{span}_vline_left_l_{left}_r_{right}.{condition}.methylation_merged_num_matrix.eps", 
        
    shell: 
        "sh scripts/plot_footprint_and_methylation.sh"
        " {input.ordered_footprint_matrix} {input.ordered_methylation_matrix}"
        " {input.site_mnase_data} {input.gnuplt_mnase_params}" 
        " {input.gnuplt_footprint_params} {output.ordered_footprint_pdf}"
        " {output.ordered_methylation_pdf}"
        " {output.ordered_footprint_gplt} {output.ordered_methylation_gplt}"
        " {wildcards.lextend} {wildcards.rextend}"
        " {wildcards.lflank} {wildcards.rflank}"
        " {input.gnuplt_methylation_params}"
        " {output.ordered_footprint_mat}"
        " {output.ordered_methylation_mat}"
        #" {output.ordered_footprint_eps}"
        #" {output.ordered_methylation_eps}"


###################### Do this for cobinding data #####################################  

rule order_extended_footprint_vec_cobinding_to_matrix:
    input:
        annotated_vec = "annotated_cobinding_status/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_sec_peak_at_{s_peak}_lex_{lextend}_rex_{rextend}_site_{site}_sam_flag_{sam_flag}_cobinding.tsv", 
        cobinding_footprint_vec_verbose = "cobinding_footprints/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_lex_{lextend}_rex_{rextend}_site_{site}_sam_flag_{sam_flag}_verbose.tsv"
    
    params:
        
    output:
        ordered_footprint_vec = "ordered_annotated_cobound_states/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_sec_peak_at_{s_peak}_lex_{lextend}_rex_{rextend}_site_{site}_sam_flag_{sam_flag}_cobinding.footprint_ordered.tsv",
        ordered_methylation_vec = "ordered_annotated_cobound_states/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_sec_peak_at_{s_peak}_lex_{lextend}_rex_{rextend}_site_{site}_sam_flag_{sam_flag}_cobinding.mvec_ordered.tsv", 
    shell:
        "sh scripts/order_methylation_and_footprint_vec.sh"
        " {input.annotated_vec} {input.cobinding_footprint_vec_verbose}" 
        " {output.ordered_footprint_vec} {output.ordered_methylation_vec}"


rule sam_flag_specific_merge_cobinding:
    input:
        ordered_99_147_footprint = "ordered_annotated_cobound_states/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_sec_peak_at_{s_peak}_lex_{lextend}_rex_{rextend}_site_{site}_sam_flag_99~147_cobinding.footprint_ordered.tsv",
        ordered_83_163_footprint = "ordered_annotated_cobound_states/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_sec_peak_at_{s_peak}_lex_{lextend}_rex_{rextend}_site_{site}_sam_flag_83~163_cobinding.footprint_ordered.tsv",
        ordered_99_147_methylation = "ordered_annotated_cobound_states/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_sec_peak_at_{s_peak}_lex_{lextend}_rex_{rextend}_site_{site}_sam_flag_99~147_cobinding.mvec_ordered.tsv", 
        ordered_83_163_methylation = "ordered_annotated_cobound_states/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_sec_peak_at_{s_peak}_lex_{lextend}_rex_{rextend}_site_{site}_sam_flag_83~163_cobinding.mvec_ordered.tsv", 
    output:
        merged_footprints = "ordered_annotated_cobound_states/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_sec_peak_at_{s_peak}_lex_{lextend}_rex_{rextend}_site_{site}_cobinding.footprint_merged_ordered.tsv", 
        merged_methylation = "ordered_annotated_cobound_states/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_sec_peak_at_{s_peak}_lex_{lextend}_rex_{rextend}_site_{site}_cobinding.mvec_merged_ordered.tsv", 
    shell:
        "sh scripts/merge_and_order_footprints_and_methylation.sh"
        " {input.ordered_99_147_footprint} {input.ordered_83_163_footprint}"
        " {input.ordered_99_147_methylation} {input.ordered_83_163_methylation}"
        " {output.merged_footprints} {output.merged_methylation}" 
        
rule plot_merged_reads_footprint_and_methylation_together_cobinding: 
    input:
        ordered_footprint_matrix = "ordered_annotated_cobound_states/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_sec_peak_at_{s_peak}_lex_{lextend}_rex_{rextend}_site_{site}_cobinding.footprint_merged_ordered.tsv",
        site_mnase_data = "site_specific_mnase_data/taken_from_{ref_flank}_cluster_{cl_id}_span_from_center_{span}_vline_left_{left}_right_{right}_site_{site}.tsv", 
        ordered_methylation_matrix = "ordered_annotated_cobound_states/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_sec_peak_at_{s_peak}_lex_{lextend}_rex_{rextend}_site_{site}_cobinding.mvec_merged_ordered.tsv", 

        gnuplt_mnase_params = "utils/gnuplot_base_files/manse_params.gplt", 
        gnuplt_footprint_params = "utils/gnuplot_base_files/footprint_params.gplt", 
        gnuplt_methylation_params = "utils/gnuplot_base_files/methylation_params.gplt", 
 
    params:
    output:
        ordered_footprint_pdf = "plots/extended_gnuplot_based_figures_cobinding/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_sec_peak_at_{s_peak}_lex_{lextend}_rex_{rextend}_site_{site}_mnase_ref_{ref_flank}_span_{span}_vline_left_l_{left}_r_{right}_cobinding.footprint_merged_num_matrix.pdf", 
        ordered_methylation_pdf = "plots/extended_gnuplot_based_figures_cobinding/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_sec_peak_at_{s_peak}_lex_{lextend}_rex_{rextend}_site_{site}_mnase_ref_{ref_flank}_span_{span}_vline_left_l_{left}_r_{right}_cobinding.methylation_merged_num_matrix.pdf", 
        ordered_footprint_gplt = "plots/extended_gnuplot_based_figures_cobinding/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_sec_peak_at_{s_peak}_lex_{lextend}_rex_{rextend}_site_{site}_mnase_ref_{ref_flank}_span_{span}_vline_left_l_{left}_r_{right}_cobinding.footprint_merged_num_matrix.gplt", 
        ordered_methylation_gplt = "plots/extended_gnuplot_based_figures_cobinding/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_sec_peak_at_{s_peak}_lex_{lextend}_rex_{rextend}_site_{site}_mnase_ref_{ref_flank}_span_{span}_vline_left_l_{left}_r_{right}_cobinding.methylation_merged_num_matrix.gplt",
        ordered_footprint_mat = "plots/extended_gnuplot_based_figures_cobinding/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_sec_peak_at_{s_peak}_lex_{lextend}_rex_{rextend}_site_{site}_mnase_ref_{ref_flank}_span_{span}_vline_left_l_{left}_r_{right}_cobinding.footprint_merged_num_matrix.mat.tsv", 
        ordered_methylation_mat = "plots/extended_gnuplot_based_figures_cobinding/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_sec_peak_at_{s_peak}_lex_{lextend}_rex_{rextend}_site_{site}_mnase_ref_{ref_flank}_span_{span}_vline_left_l_{left}_r_{right}_cobinding.methylation_merged_num_matrix.mat.tsv", 
        #ordered_footprint_eps = "plots/extended_gnuplot_based_figures/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_extend_from_peak_left_{lextend}_right_{rextend}_mnase_ref_{ref_flank}_span_{span}_vline_left_l_{left}_r_{right}.{condition}.footprint_merged_num_matrix.eps", 
        #ordered_methylation_eps = "plots/extended_gnuplot_based_figures/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_extend_from_peak_left_{lextend}_right_{rextend}_mnase_ref_{ref_flank}_span_{span}_vline_left_l_{left}_r_{right}.{condition}.methylation_merged_num_matrix.eps", 
        
    shell: 
        #"sh scripts/plot_footprint_and_methylation.sh"
        "sh scripts/plot_footprint_and_methylation_cobinding.sh"
        " {input.ordered_footprint_matrix} {input.ordered_methylation_matrix}"
        " {input.site_mnase_data} {input.gnuplt_mnase_params}" 
        " {input.gnuplt_footprint_params} {output.ordered_footprint_pdf}"
        " {output.ordered_methylation_pdf}"
        " {output.ordered_footprint_gplt} {output.ordered_methylation_gplt}"
        " {wildcards.lextend} {wildcards.rextend}"
        " {wildcards.lflank} {wildcards.rflank}"
        " {input.gnuplt_methylation_params}"
        " {output.ordered_footprint_mat}"
        " {output.ordered_methylation_mat}"
        " {wildcards.s_peak}"
        #" {output.ordered_footprint_eps}"
        #" {output.ordered_methylation_eps}"


