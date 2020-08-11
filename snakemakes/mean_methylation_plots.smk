rule mean_methylation_plots: 
    input:
        tandem_matrix = "tandem_matrix/{bigwig_annotation1}_and_{bigwig_annotation2}_to_{bed}_tandem_matrix.csv.gz"
    params:
        system_annotation = lambda wildcards: config["system_annotation"][wildcards.bigwig_annotation1] # any one of `1` or `2` is good
    output:
        methylation_cmean_png = "plots/methylation/{bigwig_annotation1}_and_{bigwig_annotation2}_to_{bed}_cmean_map.png",
        methylation_cmean_tsv = "tandem_matrix/{bigwig_annotation1}_and_{bigwig_annotation2}_to_{bed}_cmean_map.tsv"
    shell:
        "Rscript $NGS_SCRIPTS_DIR/colmeans_gz_kmeans.R"
        " {input.tandem_matrix} {output.methylation_cmean_tsv}"
        " {output.methylation_cmean_png} {params.system_annotation}"
