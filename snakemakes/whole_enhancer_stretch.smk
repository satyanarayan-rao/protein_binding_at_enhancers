rule percentage_methylation_on_whole_enh:
    input:
        footprint_bed = "fragments_mapped_to_enhancer_center/{sample}_to_{bed}_footprint.bed", 
        
    params:
        
    output:
        whole_enh_per_methylation = "whole_enhancer_percentage_methylation/{sample}_to_{bed}_occl_and_per_methylation.tsv",
        whole_enh_per_methylation_pkl = "whole_enhancer_percentage_methylation/{sample}_to_{bed}_occl_and_per_methylation.pkl",
    shell:
        "python scripts/occluded_edges_and_percentage_methylation.py"
        " {input.footprint_bed}"
        " {output.whole_enh_per_methylation} {output.whole_enh_per_methylation_pkl}" 
