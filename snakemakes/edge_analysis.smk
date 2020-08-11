rule occluded_edges_on_reads:
    input:
        footprint_capped_bed = "footprint_min_size_selection/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_tf_footprint_min_length_10_wobble_gap_1.bed"
    params:
    output:
        occluded_bed = "occluded_edges_on_methylation_vec/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_min_flen_10_wobble_gap_1_occluded.tsv",
        occluded_pkl = "occluded_edges_on_methylation_vec/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_min_flen_10_wobble_gap_1_occluded.pkl",
    shell:
        "python scripts/occluded_edges_and_percentage_methylation.py"
        " {input.footprint_capped_bed}"
        " {output.occluded_bed} {output.occluded_pkl}" 

