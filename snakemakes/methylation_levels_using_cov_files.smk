rule methylation_report_to_dict: 
    input:
        report_file = "coverage_to_cytosine/{sample}_report.txt.gz"
    params:
    output:
        pkl_file = "coverage_or_report_dict/{sample}_report.pkl"
    shell:
        "python scripts/report2pkl.py {input.report_file} {output.pkl_file}"

rule generate_methylation_level_matrix_using_cov_files: 
    input:
        cpg_cov_report_file = lambda wildcards: config["report_dict_annotation"][wildcards.report_annotation_cpg],
        gpc_cov_report_file = lambda wildcards: config["report_dict_annotation"][wildcards.report_annotation_gpc],
        input_bed = "input_bed/{bed}.bed"
        #fasta_out = "fasta_files/reordered_20clust_147_50PNE_open_v1_peak_{flank}.fa"
    params:
    output:
        methylation_level_matrix = "methylation_level_matrix/{report_annotation_cpg}_and_{report_annotation_gpc}_for_{bed}_matrix.tsv.gz"
    shell:
        "python scripts/build_methylation_level_matrix.py {input.input_bed} {input.cpg_cov_report_file} {input.gpc_cov_report_file} {output.methylation_level_matrix}"

rule mean_methylation_plots_per_cluster: 
    input:
        methylation_level_matrix = "methylation_level_matrix/{report_annotation_cpg}_and_{report_annotation_gpc}_for_{bed}_matrix.tsv.gz"
    params:
        system_annotation = lambda wildcards: config["system_annotation"][wildcards.report_annotation_cpg], # any one of `cpg` or `gpc` is good
        figure_height = 4, 
        rolling_window_size = 50

    output:
        methylation_cmean_per_cluster_png = "plots/methylation/{report_annotation_cpg}_and_{report_annotation_gpc}_for_{bed}_mean_methylation_clust_{cl_id}_plot.png",
        frac_methyl_sites_per_cluster_png = "plots/methylation/{report_annotation_cpg}_and_{report_annotation_gpc}_for_{bed}_frac_clust_{cl_id}_plot.png",
        count_methyl_sites_per_cluster_png = "plots/methylation/{report_annotation_cpg}_and_{report_annotation_gpc}_for_{bed}_count_clust_{cl_id}_plot.png",
        methylation_cmean_per_cluster_tsv = "colmeans_of_cluster_methylation/{report_annotation_cpg}_and_{report_annotation_gpc}_for_{bed}_mean_methylation_clust_{cl_id}_plot.tsv",
        frac_methyl_sites_per_cluster_tsv = "colmeans_of_cluster_methylation/{report_annotation_cpg}_and_{report_annotation_gpc}_for_{bed}_frac_clust_{cl_id}_plot.tsv",
        count_methyl_sites_per_cluster_tsv = "colmeans_of_cluster_methylation/{report_annotation_cpg}_and_{report_annotation_gpc}_for_{bed}_count_clust_{cl_id}_plot.tsv"
    shell:
        "Rscript $NGS_SCRIPTS_DIR/colmeans_gz_kmeans_one_cluster.R"
        " {input.methylation_level_matrix}"
        " {wildcards.cl_id}"
        " {output.methylation_cmean_per_cluster_png}"
        " {output.frac_methyl_sites_per_cluster_png}"
        " {output.count_methyl_sites_per_cluster_png}"
        " {output.methylation_cmean_per_cluster_tsv}"
        " {output.frac_methyl_sites_per_cluster_tsv}"
        " {output.count_methyl_sites_per_cluster_tsv}"
        " {params.system_annotation} {params.figure_height}"
        " {params.rolling_window_size}"
