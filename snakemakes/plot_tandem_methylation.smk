rule plot_tandem_methylation: 
    input:
        tandem_matrix = "tandem_matrix/{bigwig_annotation1}_and_{bigwig_annotation2}_to_{bed}_tandem_matrix.csv.gz",
        base_gnuplt_methylation = "utils/gnuplot_base_files/methylation_plot.gplt"
    params:
    output:
        methylation_heatmap_eps = "plots/methylation/{bigwig_annotation1}_and_{bigwig_annotation2}_to_{bed}_tandem_matrix.eps", 
        methylation_heatmap_pdf = "plots/methylation/{bigwig_annotation1}_and_{bigwig_annotation2}_to_{bed}_tandem_matrix.pdf",
        methylation_heatmap_plt = "plots/methylation/{bigwig_annotation1}_and_{bigwig_annotation2}_to_{bed}_tandem_matrix.gplt"
    
    shell:
        "sh scripts/plot_methylation.sh"
        " {input.tandem_matrix} {input.base_gnuplt_methylation}"
        " {output.methylation_heatmap_eps} {output.methylation_heatmap_pdf}"
        " {output.methylation_heatmap_plt}"
