rule plot_percentage_methylation:  # `a`: open enhancers, `b`: close enhancers 
    input:
        occluded_tsv_a = "occluded_edges_on_methylation_vec/{sample}_to_{bed_a}_lf_{lflank}_rf_{rflank}_min_flen_10_wobble_gap_1_occluded.tsv",
        occluded_tsv_b = "occluded_edges_on_methylation_vec/{sample}_to_{bed_b}_lf_{lflank}_rf_{rflank}_min_flen_10_wobble_gap_1_occluded.tsv"
    params:
        percentage_field = 4, 
        bw = 5
    output:
        dist_a_vs_b_pdf = "plots/percentage_methyation_on_reads/{sample}_to_{bed_a}_vs_{bed_b}_lf_{lflank}_rf_{rflank}.dist.pdf", 
        dist_a_vs_b_png = "plots/percentage_methyation_on_reads/{sample}_to_{bed_a}_vs_{bed_b}_lf_{lflank}_rf_{rflank}.dist.png",
        smooth_dist_a_vs_b_pdf = "plots/percentage_methyation_on_reads/{sample}_to_{bed_a}_vs_{bed_b}_lf_{lflank}_rf_{rflank}.dist_smooth.pdf", 
        smooth_dist_a_vs_b_png = "plots/percentage_methyation_on_reads/{sample}_to_{bed_a}_vs_{bed_b}_lf_{lflank}_rf_{rflank}.dist_smooth.png",
        smooth_dist_a_vs_b_ecdf_pdf = "plots/percentage_methyation_on_reads/{sample}_to_{bed_a}_vs_{bed_b}_lf_{lflank}_rf_{rflank}.ecdf_smooth.pdf", 
        smooth_dist_a_vs_b_ecdf_png = "plots/percentage_methyation_on_reads/{sample}_to_{bed_a}_vs_{bed_b}_lf_{lflank}_rf_{rflank}.ecdf_smooth.png", 
        smooth_dist_a_vs_b_cnt_pdf = "plots/percentage_methyation_on_reads/{sample}_to_{bed_a}_vs_{bed_b}_lf_{lflank}_rf_{rflank}.cnt_smooth.pdf", 
        smooth_dist_a_vs_b_cnt_png = "plots/percentage_methyation_on_reads/{sample}_to_{bed_a}_vs_{bed_b}_lf_{lflank}_rf_{rflank}.cnt_smooth.png" 
        
    shell:
        "Rscript scripts/plot_percentage_methylation_on_reads.R"
        " {input.occluded_tsv_a} {input.occluded_tsv_b}"
        " {params.percentage_field} {params.bw}"
        " {output.dist_a_vs_b_pdf} {output.dist_a_vs_b_png}"  
        " {output.smooth_dist_a_vs_b_pdf} {output.smooth_dist_a_vs_b_png}"
        " {output.smooth_dist_a_vs_b_ecdf_pdf} {output.smooth_dist_a_vs_b_ecdf_png}"  
        " {output.smooth_dist_a_vs_b_cnt_pdf} {output.smooth_dist_a_vs_b_cnt_png}"

rule plot_percentage_methylation_for_whole_enh:  # `a`: open enhancers, `b`: close enhancers 
    input:
        occluded_tsv_a = "whole_enhancer_percentage_methylation/{sample}_to_{bed_a}_occl_and_per_methylation.tsv",
        occluded_tsv_b = "whole_enhancer_percentage_methylation/{sample}_to_{bed_b}_occl_and_per_methylation.tsv"
    params:
        percentage_field = 4, 
        bw = 5
    output:
        dist_a_vs_b_pdf = "plots/whole_enh_percentage_methyation_on_reads/{sample}_to_{bed_a}_vs_{bed_b}.dist.pdf", 
        dist_a_vs_b_png = "plots/whole_enh_percentage_methyation_on_reads/{sample}_to_{bed_a}_vs_{bed_b}.dist.png",
        smooth_dist_a_vs_b_pdf = "plots/whole_enh_percentage_methyation_on_reads/{sample}_to_{bed_a}_vs_{bed_b}.dist_smooth.pdf", 
        smooth_dist_a_vs_b_png = "plots/whole_enh_percentage_methyation_on_reads/{sample}_to_{bed_a}_vs_{bed_b}.dist_smooth.png",
        smooth_dist_a_vs_b_ecdf_pdf = "plots/whole_enh_percentage_methyation_on_reads/{sample}_to_{bed_a}_vs_{bed_b}.ecdf_smooth.pdf", 
        smooth_dist_a_vs_b_ecdf_png = "plots/whole_enh_percentage_methyation_on_reads/{sample}_to_{bed_a}_vs_{bed_b}.ecdf_smooth.png", 
        smooth_dist_a_vs_b_cnt_pdf = "plots/whole_enh_percentage_methyation_on_reads/{sample}_to_{bed_a}_vs_{bed_b}.cnt_smooth.pdf", 
        smooth_dist_a_vs_b_cnt_png = "plots/whole_enh_percentage_methyation_on_reads/{sample}_to_{bed_a}_vs_{bed_b}.cnt_smooth.png" 
        
    shell:
        "Rscript scripts/plot_percentage_methylation_on_reads.R"
        " {input.occluded_tsv_a} {input.occluded_tsv_b}"
        " {params.percentage_field} {params.bw}"
        " {output.dist_a_vs_b_pdf} {output.dist_a_vs_b_png}"  
        " {output.smooth_dist_a_vs_b_pdf} {output.smooth_dist_a_vs_b_png}"
        " {output.smooth_dist_a_vs_b_ecdf_pdf} {output.smooth_dist_a_vs_b_ecdf_png}"  
        " {output.smooth_dist_a_vs_b_cnt_pdf} {output.smooth_dist_a_vs_b_cnt_png}"

