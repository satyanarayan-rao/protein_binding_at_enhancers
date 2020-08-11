rule map_5mc_lvels_to_enhancers: 
    input:
        bed_file = "input_bed/{bed}.bed",
        bw_file = lambda wildcards: config["bigwig_annotation"][wildcards.bigwig_annotation]
    params:
        strand_column = 6
    output:
        raw_signal = "bw_map/{bigwig_annotation}_to_{bed}_map.csv.gz", 
        eom = "bw_map/{bigwig_annotation}_to_{bed}_map.csv.gz_e.csv.gz", # eom: enrichment over mean
        exact_signal = "bw_map/{bigwig_annotation}_to_{bed}_map.csv.gz_exact_stats.tsv"
        
    shell:
        "python $NGS_SCRIPTS_DIR/map_bw_to_bed_strand_aware.py"
        " {input.bw_file} {input.bed_file} {output.raw_signal} {params.strand_column}"
