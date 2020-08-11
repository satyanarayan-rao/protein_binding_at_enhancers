rule map_mnase: 
    input: 
        bed_file = "input_bed/{bed}.bed", 
        bw_file = "bigwigs/{mnase_signal}.bigwig"
    params:
        strand_column = 6
    output:
        raw_signal = "mnase_map/{mnase_signal}_to_{bed}_map.csv.gz", 
        eom = "mnase_map/{mnase_signal}_to_{bed}_map.csv.gz_e.csv.gz",
        exact_signal = "mnase_map/{mnase_signal}_to_{bed}_map.csv.gz_exact_stats.tsv"
    shell:
        "python $NGS_SCRIPTS_DIR/map_bw_to_bed_strand_aware.py"
        " {input.bw_file} {input.bed_file} {output.raw_signal} {params.strand_column}"
rule mnase_heatmap:
    input: 
        raw_signal = "mnase_map/{mnase_signal}_to_{bed}_map.csv.gz", 
        base_gnuplt_raw_signal = "utils/gnuplot_base_files/mnase_raw.gplt", 
        eom = "mnase_map/{mnase_signal}_to_{bed}_map.csv.gz_e.csv.gz",
        base_gnuplt_eom_signal = "utils/gnuplot_base_files/mnase_eom.gplt",
        exact_signal = "mnase_map/{mnase_signal}_to_{bed}_map.csv.gz_exact_stats.tsv"
    params: 
        flank_str = lambda wildcards: config["bed_flank_annotation"][wildcards.bed]
    output:
        raw_eps = "plots/mnase/{mnase_signal}_to_{bed}_raw.eps",
        raw_pdf = "plots/mnase/{mnase_signal}_to_{bed}_raw.pdf",
        raw_gplt = "plots/mnase/{mnase_signal}_to_{bed}_raw.gplt", 
        eom_eps = "plots/mnase/{mnase_signal}_to_{bed}_eom.eps",
        eom_pdf = "plots/mnase/{mnase_signal}_to_{bed}_eom.pdf",
        eom_gplt = "plots/mnase/{mnase_signal}_to_{bed}_eom.gplt"
    shell:
        "sh scripts/plot_mnase.sh {input.raw_signal} {input.base_gnuplt_raw_signal}"
        " {input.eom} {input.base_gnuplt_eom_signal} {input.exact_signal}"
        " {output.raw_eps} {output.raw_gplt}"
        " {output.eom_eps} {output.eom_gplt}" 
        " {params.flank_str}"
        " {output.raw_pdf} {output.eom_pdf}"
rule mean_mnase_plots:
    input:
        raw_signal = "mnase_map/{mnase_signal}_to_{bed}_map.csv.gz", 
        eom_signal = "mnase_map/{mnase_signal}_to_{bed}_map.csv.gz_e.csv.gz" 
    params:
        rollmean_window = 50 
    output:
        raw_rolled_cmean_png = "plots/mnase/{mnase_signal}_to_{bed}_raw_rolled_cmean_map.png", 
        eom_rolled_cmean_png = "plots/mnase/{mnase_signal}_to_{bed}_eom_rolled_cmean_map.png", 
        raw_cmean_tsv = "mnase_map/{mnase_signal}_to_{bed}_raw_cmean.tsv", 
        raw_rolled_cmean_tsv = "mnase_map/{mnase_signal}_to_{bed}_raw_rolled_cmean.tsv", 
        eom_cmean_tsv = "mnase_map/{mnase_signal}_to_{bed}_eom_cmean.tsv", 
        eom_rolled_cmean_tsv = "mnase_map/{mnase_signal}_to_{bed}_eom_rolled_cmean.tsv", 
    shell:        
        "Rscript $NGS_SCRIPTS_DIR/colmeans_plus_rolling_mean_gz_kmeans.R"
        " {input.raw_signal} {output.raw_cmean_tsv} {params.rollmean_window} {output.raw_rolled_cmean_tsv}"
        " {output.raw_rolled_cmean_png} \"{wildcards.mnase_signal}\""
        "; Rscript $NGS_SCRIPTS_DIR/colmeans_plus_rolling_mean_gz_kmeans.R"
        " {input.eom_signal} {output.eom_cmean_tsv} {params.rollmean_window} {output.eom_rolled_cmean_tsv}"
        " {output.eom_rolled_cmean_png} \"{wildcards.mnase_signal}\""

rule mean_mnase_plots_per_cluster:
    input:
        raw_signal = "mnase_map/{mnase_signal}_to_{bed}_map.csv.gz", 
        eom_signal = "mnase_map/{mnase_signal}_to_{bed}_map.csv.gz_e.csv.gz" 
    params:
        rollmean_window = 50, 
        figure_height = 4
    output:
        raw_rolled_cmean_png = "plots/mnase/{mnase_signal}_to_{bed}_raw_rolled_cmean_clust_{cl_id}_map.png", 
        raw_cmean_and_rolled_tsv = "mnase_map/{mnase_signal}_to_{bed}_raw_cmean_clust_{cl_id}.tsv", 
        eom_rolled_cmean_png = "plots/mnase/{mnase_signal}_to_{bed}_eom_rolled_cmean_clust_{cl_id}_map.png", 
        eom_cmean_and_rolled_tsv = "mnase_map/{mnase_signal}_to_{bed}_eom_cmean_clust_{cl_id}.tsv", 
    shell:        
        "Rscript $NGS_SCRIPTS_DIR/colmeans_plus_rolling_mean_gz_kmeans_for_cluster.R"
        " {input.raw_signal} {wildcards.cl_id}"
        " {output.raw_rolled_cmean_png}"
        " {output.raw_cmean_and_rolled_tsv}"
        "  \"{wildcards.mnase_signal}\""
        " {params.figure_height} {params.rollmean_window}"
        " \"Mean MNase (raw)\""
        "; Rscript $NGS_SCRIPTS_DIR/colmeans_plus_rolling_mean_gz_kmeans_for_cluster.R"
        " {input.eom_signal} {wildcards.cl_id}"
        " {output.eom_rolled_cmean_png}"
        " {output.eom_cmean_and_rolled_tsv}"
        "  \"{wildcards.mnase_signal}\""
        " {params.figure_height} {params.rollmean_window}"
        " \"Mean MNase (e.o.m)\""

rule peak_locations_per_site: 
    input:
        eom_signal = "mnase_map/{mnase_signal}_to_{bed}_map.csv.gz_e.csv.gz" 
    params:
        rollmean_window = 50
    output:
        site_wise_eom_plot = "site_wise_plots/{mnase_signal}_to_{bed}_site_{site}_eom.pdf", 
        site_wise_eom_rmean_plot = "site_wise_plots/{mnase_signal}_to_{bed}_site_{site}_eom_rmean.pdf",
        site_wise_max_peak_locations = "site_wise_plots/{mnase_signal}_to_{bed}_site_{site}_peak_location.tsv"
    shell:
        "Rscript scripts/site_wise_mnase_signal.R {input.eom_signal}"
        " \"{wildcards.site}\" {output.site_wise_eom_plot}"
        " {output.site_wise_eom_rmean_plot}"
        " {output.site_wise_max_peak_locations}"  
        " {params.rollmean_window}"
