rule get_enhancer_site_sequences:
    input: 
        input_bed = "input_bed/reordered_20clust_147_50PNE_open_v1_peak_{flank}.bed",
        genome_fa = "~/data/ucsc/dm/dm3/dm3.fa"
    params:
    output:
        fasta_out = "fasta_files/reordered_20clust_147_50PNE_open_v1_peak_{flank}.fa"
    shell:
        "sh scripts/get_fasta.sh {input.input_bed} {input.genome_fa} {output.fasta_out}"

