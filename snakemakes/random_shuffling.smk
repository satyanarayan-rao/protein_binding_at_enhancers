rule get_methylated_strings: 
    input:
        bam_file = "bismark_mapped/{sample}.bam"
    params:
    output:
        read_level_methylation_status = "bisulphite_sequences/{sample}_bisulphite.bed.gz"
    shell:
        "sh scripts/get_methylation_sequences.sh {input.bam_file} {output.read_level_methylation_status}"
rule select_overlapping_or_adjacent: 
    input:
        read_level_methylation_status = "bisulphite_sequences/{sample}_bisulphite.bed.gz"
    params:
    output:
        overlapping_or_adjacent = "bisulphite_overlapping_or_adjacent/{sample}_overlapping_or_adjacent.bisulphite.bed.gz"
    shell:
        "zcat {input.read_level_methylation_status}"
        " | egrep \"_overlapping|_adjacent\" "
        " | gzip - > {output.overlapping_or_adjacent}"


rule get_genomic_seq_for_overlap_or_adj:
    input:
        overlapping_or_adjacent = "bisulphite_overlapping_or_adjacent/{sample}_overlapping_or_adjacent.bisulphite.bed.gz"
    params:
    output:
        genomic_seq_zero_based = "bisulphite_overlapping_or_adjacent/{sample}_overlapping_or_adjacent.genomic.0-based.bed.gz",
        genomic_seq_one_based = "bisulphite_overlapping_or_adjacent/{sample}_overlapping_or_adjacent.genomic.1-based.bed.gz"
    shell:
        "sh scripts/get_genomic_seq.sh {input.overlapping_or_adjacent} {output.genomic_seq_zero_based} {output.genomic_seq_one_based}"


rule bisulphite_seq_mapped_to_enhancers: 
    input:
        overlapping_or_adjacent = "bisulphite_overlapping_or_adjacent/{sample}_overlapping_or_adjacent.bisulphite.bed.gz",
        target_bed = "input_bed/{bed}.bed"
    params:
        chunk_size = 1000000 
    output:
        mapped_to_enhancer = "fragments_mapped_to_enhancer_center/{sample}_to_{bed}_mapped.bsseq.bed" # bsseq : BiSuphite Sequence
    shell:
        "Rscript scripts/intersection.R {input.target_bed}"
        " {input.overlapping_or_adjacent} {params.chunk_size}"
        " {output.mapped_to_enhancer}"  
        
rule genomic_seq_mapped_to_enhancers: 
    input:
        mapped_to_enhancer = "fragments_mapped_to_enhancer_center/{sample}_to_{bed}_mapped.bsseq.bed" # bsseq : BiSuphite Sequence
    output:
        genomic_seq = "fragments_mapped_to_enhancer_center/{sample}_to_{bed}_mapped.genomic.1-based.bed" # genomic: Genomic Sequence
    shell:
        "sh scripts/gemomic_seq_for_bsseq.sh {input.mapped_to_enhancer} {output.genomic_seq} "
  
rule random_shuffling:
    input:
        genomic_seq = "fragments_mapped_to_enhancer_center/{sample}_to_{bed}_mapped.genomic.1-based.bed", # genomic: Genomic Sequence
        mapped_to_enhancer_mvec = "fragments_mapped_to_enhancer_center/{sample}_to_{bed}_mapped.bed"
    params:
    output:
        random_mvec = "fragments_mapped_to_enhancer_center_shuffled/{sample}_to_{bed}_seed_{seed}_mapped.bed", 
        random_mvec_verbose = "fragments_mapped_to_enhancer_center_shuffled/{sample}_to_{bed}_seed_{seed}_mapped.verbose.bed", 
    shell:
        "python scripts/shuffle_vec.py {input.genomic_seq} {input.mapped_to_enhancer_mvec} {output.random_mvec} {output.random_mvec_verbose} {wildcards.seed}" 

rule footprints_on_shuffled_vec: 
    input:
        random_mvec = "fragments_mapped_to_enhancer_center_shuffled/{sample}_to_{bed}_seed_{seed}_mapped.bed"
    params:
    output:
        footprints_on_shuf = "footprints_on_shuffled_mvec/{sample}_to_{bed}_seed_{seed}_footprint.bed", 
        footprints_cnt_matrix = "footprints_on_shuffled_mvec/{sample}_to_{bed}_seed_{seed}_footprint.cnt", 
    shell:
        "python scripts/call_footprints_v2.py"
        " {input.random_mvec} {output.footprints_on_shuf} {output.footprints_cnt_matrix}"   
rule footprints_on_observed:
    input:
        observed = "fragments_mapped_to_enhancer_center/{sample}_to_{bed}_mapped.bed"
    params:
    output:
        obs_footprints = "observed_footprints/{sample}_to_{bed}_footprint.bed",  
        obs_footprints_cnt_matrix = "observed_footprints/{sample}_to_{bed}_footprint.cnt",  
    shell:
        "python scripts/call_footprints_v2.py"
        " {input.observed} {output.obs_footprints} {output.obs_footprints_cnt_matrix}"

def get_shuf_files_open_enh (wildcards):
    flist = [] 
    for shuf in config["random_shuffle_seed"]["seeds"]:
        flist.append("footprints_on_shuffled_mvec/merged_S2_to_reordered_20clust_147_50PNE_open_v1_peak_seed_%s_footprint.cnt"%(shuf)) # append for open enhancers
    to_return = flist
    return to_return
def get_shuf_files_closed_enh (wildcards):
    flist = [] 
    for shuf in config["random_shuffle_seed"]["seeds"]:
        flist.append("footprints_on_shuffled_mvec/merged_S2_to_50PNE_closed_peak_seed_%s_footprint.cnt"%(shuf)) # append for open enhancers
    to_return = flist
    return to_return
rule prepare_files_for_hist: 
    input:
        open_enh = "observed_footprints/merged_S2_to_reordered_20clust_147_50PNE_open_v1_peak_footprint.cnt",
        closed_enh = "observed_footprints/merged_S2_to_50PNE_closed_peak_footprint.cnt",
        shuf_files_open_enh =  lambda wildcards: get_shuf_files_open_enh (wildcards),
        shuf_files_closed_enh = lambda wildcards: get_shuf_files_closed_enh (wildcards)
    params:
        rseed_str = "@".join(map(str,config["random_shuffle_seed"]["seeds"]))
    output:
        hist_file = "observed_vs_random_footprint/observed_vs_random.hist.tsv" 
    shell:
        "sh scripts/prepare_hist_file.sh {input.open_enh} {input.closed_enh}"
        " \"{input.shuf_files_open_enh}\" \"{input.shuf_files_closed_enh}\" {params.rseed_str}"
        " {output.hist_file}" 
def get_shuf_files_open_enh_whole (wildcards):
    flist = [] 
    for shuf in config["random_shuffle_seed"]["seeds"]:
        flist.append("footprints_on_shuffled_mvec/merged_S2_to_open_S2_STARRSEQ_formatted_seed_%s_footprint.cnt"%(shuf)) # append for open enhancers
    to_return = flist
    return to_return
def get_shuf_files_closed_enh_whole (wildcards):
    flist = [] 
    for shuf in config["random_shuffle_seed"]["seeds"]:
        flist.append("footprints_on_shuffled_mvec/merged_S2_to_closed_S2_STARRSEQ_formatted_seed_%s_footprint.cnt"%(shuf)) # append for open enhancers
    to_return = flist
    return to_return
rule prepare_files_for_hist_for_whole_enh: 
    input:
        open_enh = "observed_footprints/merged_S2_to_open_S2_STARRSEQ_formatted_footprint.cnt",
        closed_enh = "observed_footprints/merged_S2_to_closed_S2_STARRSEQ_formatted_footprint.cnt",
        shuf_files_open_enh =  lambda wildcards: get_shuf_files_open_enh_whole (wildcards),
        shuf_files_closed_enh = lambda wildcards: get_shuf_files_closed_enh_whole (wildcards)
    params:
        rseed_str = "@".join(map(str,config["random_shuffle_seed"]["seeds"]))
    output:
        hist_file = "observed_vs_random_footprint/observed_vs_random_whole_enh.hist.tsv" 
    shell:
        "sh scripts/prepare_hist_file.sh {input.open_enh} {input.closed_enh}"
        " \"{input.shuf_files_open_enh}\" \"{input.shuf_files_closed_enh}\" {params.rseed_str}"
        " {output.hist_file}"
