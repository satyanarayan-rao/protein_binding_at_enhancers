rule build_tandem_matrix: 
    input: 
        cpg_matrix = lambda wildcards: "bw_map/{bw}_to_{bed}_map.csv.gz".format(bw = wildcards.bigwig_annotation1, bed = wildcards.bed), 
        gpc_matrix = lambda wildcards: "bw_map/{bw}_to_{bed}_map.csv.gz".format(bw = wildcards.bigwig_annotation2, bed = wildcards.bed)
    params:
    output:
        tandem_matrix = "tandem_matrix/{bigwig_annotation1}_and_{bigwig_annotation2}_to_{bed}_tandem_matrix.csv.gz"
    shell:
        "Rscript scripts/add_two_gzip_matrix.R {input.cpg_matrix} {input.gpc_matrix} {output.tandem_matrix}"
