import pandas as pd
rule intersect_bam_and_bed: 
    input:
        bam_file = "bismark_mapped/{sample}_pe_sorted.bam", 
        bed_file = "input_bed/{bed}.bed"
    params:
    output:
        intersect_file = "reads_mapped_to_enhancers/{sample}_to_{bed}_intersect.bed"
    shell:
        "bedtools intersect -abam {input.bam_file} -b {input.bed_file} -wa -wb -bed > {output.intersect_file}"
rule select_reads_from_bam: 
    input: 
        intersect_file = "reads_mapped_to_enhancers/{sample}_to_{bed}_intersect.bed",
        bam_file = "bismark_mapped/{sample}_pe.bam"
    params:
        picard_jar = "/cluster/software/modules-sw/picard-tools/2.20.1/picard.jar"
    output:
        selected_reads = "bam_subset_for_enhancers/{sample}_to_{bed}_intersect.bam"
    shell: 
        "sh scripts/picard_select_reads.sh"
        " {input.intersect_file} {input.bam_file} {output.selected_reads}"
        " {wildcards.sample} {wildcards.bed} {params.picard_jar}"

rule create_read_methylation_status_dict: 
    input:
        selected_reads = "bam_subset_for_enhancers/{sample}_to_{bed}_intersect.bam"
    params:
    output: 
        methylation_dict = "methylation_status_at_read_bases/{sample}_to_{bed}_methylation_dict.pkl"
    shell: 
        "sh scripts/create_methylation_level_dict.sh {input.selected_reads} {output.methylation_dict} {wildcards.sample}" 

rule populate_methylation_status_near_peak: 
    input: 
        intersect_file = "reads_mapped_to_enhancers/{sample}_to_{bed}_intersect.bed",
        methylation_dict = "methylation_status_at_read_bases/{sample}_to_{bed}_methylation_dict.pkl",
    params:
    output:
        methylation_status_matrix = "methylation_status_matrix_at_enhancers/{sample}_to_{bed}_slop_{slop}.tsv"
    shell:
        "python scripts/build_methylation_matrix_for_enhancers.py {input.intersect_file} {input.methylation_dict} {output.methylation_status_matrix} {wildcards.slop}"   

rule select_exact_padding: 
    input:
        methylation_status_matrix = "methylation_status_matrix_at_enhancers/{sample}_to_{bed}_slop_{slop}.tsv"
    params:
    output:
        subset_matrix = "exact_padding_from_peak_center/{sample}_to_{bed}_slop_{slop}_subset_pad_{pad}.tsv"
    shell: 
        "sh scripts/select_for_padding.sh"
        " {input.methylation_status_matrix} {output.subset_matrix} {wildcards.pad}"

rule bam2fragment_level_methylation_bedgz: 
    input:
        bam_file = "bismark_mapped/{sample}.bam"
    params:
    output:
        read_level_methylation_status = "fragment_level_methylation/{sample}_fragment_methylation_vec.bed.gz"
    shell:
        "sh scripts/get_methylation_vec.sh {input.bam_file} {output.read_level_methylation_status}"
    
rule select_for_overlapping_or_adjacent_fragments: 
    input:
        read_level_methylation_status = "fragment_level_methylation/{sample}_fragment_methylation_vec.bed.gz"
    params:
    output:
        overlapping_or_adjacent = "overlapping_or_adjacent_fragments/{sample}_overlapping_or_adjacent.bed.gz"
    shell:
        "zcat {input.read_level_methylation_status}"
        " | egrep \"_overlapping|_adjacent\" "
        " | gzip - > {output.overlapping_or_adjacent}"

rule map_fragment_level_methylation_bedgz_to_bed: 
    input:
        overlapping_or_adjacent = "overlapping_or_adjacent_fragments/{sample}_overlapping_or_adjacent.bed.gz",
        target_bed = "input_bed/{bed}.bed"
    params:
        chunk_size = 1000000 
    output:
        mapped_to_enhancer = "fragments_mapped_to_enhancer_center/{sample}_to_{bed}_mapped.bed" 
    shell:
        "Rscript scripts/intersection.R {input.target_bed}"
        " {input.overlapping_or_adjacent} {params.chunk_size}"
        " {output.mapped_to_enhancer}"
#rule assign_sam_flag: 
#    input:
#        mapped_to_enhancer = "fragments_mapped_to_enhancer_center/{sample}_to_{bed}_mapped.bed", 
#        samflag_dict = "bismark_mapped/alignment_samflag_{sample}.pkl"
#    params:
#    output:
#        mapped_with_sam_flag = "fragments_mapped_to_enhancer_center/{sample}_to_{bed}_mapped.sam.bed"
#    shell: 
#        "python scripts/assign_sam_flag.py {input.mapped_to_enhancer}"
#        " {input.samflag_dict} {output.mapped_with_sam_flag}"

rule call_footprint_on_complete_data: 
    input:
        mapped_to_enhancer = "fragments_mapped_to_enhancer_center/{sample}_to_{bed}_mapped.bed" 
    output:
        footprint_bed = "fragments_mapped_to_enhancer_center/{sample}_to_{bed}_footprint.bed", 
        footprint_cnt = "fragments_mapped_to_enhancer_center/{sample}_to_{bed}_footprint.cnt"
    shell:
        "python scripts/call_footprints_v2_no_buffer.py {input.mapped_to_enhancer}"
        " {output.footprint_bed} {output.footprint_cnt}"  
rule generate_footprint_dict:
    input:
        footprint_bed = "fragments_mapped_to_enhancer_center/{sample}_to_{bed}_footprint.bed", 
    params:
    output:
        footprint_pkl = "fragments_mapped_to_enhancer_center/{sample}_to_{bed}_footprint.pkl"
    shell:
        "python scripts/create_footprint_dict.py {input.footprint_bed} {output.footprint_pkl}" 



rule select_for_flank_from_peak_center:
    input:
        mapped_to_enhancer = "fragments_mapped_to_enhancer_center/{sample}_to_{bed}_mapped.bed" 
    params:
    output:
        selected_flank = "flank_with_methyation_status/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_mapped.bed"
    shell:
        "sh scripts/select_for_flank.sh {input.mapped_to_enhancer}"
        " {wildcards.cl_id} {wildcards.lflank} {wildcards.rflank}"
        " {output.selected_flank}"
rule select_for_flank_from_peak_center_no_cluster:
    input:
        mapped_to_enhancer = "fragments_mapped_to_enhancer_center/{sample}_to_{bed}_mapped.bed" 
    params:
    output:
        selected_flank = "flank_with_methyation_status/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_mapped.bed"
    shell:
        #"sh scripts/select_for_flank_no_cluster.sh {input.mapped_to_enhancer}" # doesn't consider strand direction
        "python scripts/select_for_flank_no_cluster.py {input.mapped_to_enhancer}" #  consider strand direction
        " {wildcards.lflank} {wildcards.rflank}"
        " {output.selected_flank}" 
#rule define_methylation_matrix: 
#    input:
#        selected_flank = "flank_with_methyation_status/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_mapped.bed",
#        samflag_dict = "bismark_mapped/alignment_samflag_{sample}.pkl"
#    params:
#    output:
#        methylation_matrix = "flank_methylation_matrix/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_matrix.bed"
#    shell:
#        "python scripts/build_flank_methylation_matrix.py {input.selected_flank}"
#        " {input.samflag_dict} {wildcards.lflank} {wildcards.rflank}"
#        " {output.methylation_matrix}" 

rule define_footprints:
    input:
        #selected_flank = "flank_with_methyation_status/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_mapped.bed",
        selected_flank = "flank_with_methyation_status/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_mapped.bed",  # keeping this agnostic to cluster id wildcard (see prevoious rule) because I want to keep it general
        #methylation_matrix = "flank_methylation_matrix/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_matrix.bed"
    params:
    output:
        #tf_footprint_on_reads = "footprints/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_tf_footprint.bed",
        #footprint_cnt_matrix = "footprints/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_tf_footprint.count"
        tf_footprint_on_reads = "footprints/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_tf_footprint.bed",
        footprint_cnt_matrix = "footprints/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_tf_footprint.count",

    shell:
        #"python scripts/find_footprints_in_methylated_reads.py"
        #"python scripts/call_footprints_v2.py"
        "python scripts/call_footprints_v2_no_buffer.py"
        " {input.selected_flank} {output.tf_footprint_on_reads} {output.footprint_cnt_matrix}"   

rule fix_wobble_footprint:
    input:
        tf_footprint_on_reads = "footprints/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_tf_footprint.bed",
    output:
        wobble_fixed_on_reads = "wobble_fixed_footprints/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_wobble_fixed_footprint_gap_{gap}.bed"
    shell:
        "python scripts/fix_wobble.py {input.tf_footprint_on_reads}"
        " {wildcards.gap} {output.wobble_fixed_on_reads}" 
        
rule min_cap_footprint_size: 
    input: 
        wobble_fixed_on_reads = "wobble_fixed_footprints/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_wobble_fixed_footprint_gap_{gap}.bed"
    output:
        footprint_capped = "footprint_min_size_selection/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_tf_footprint_min_length_{min_len}_wobble_gap_{gap}.bed"
    shell:
        "python scripts/keep_all_footprint_ge_th.py {input.wobble_fixed_on_reads}"
        " {wildcards.min_len} {output.footprint_capped}"

rule generate_min_cap_footprint_dict:
    input:
        footprint_capped = "footprint_min_size_selection/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_tf_footprint_min_length_{min_len}_wobble_gap_{gap}.bed"
    params:
    output:
        footprint_capped_pkl = "footprint_min_size_selection/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_tf_footprint_min_length_{min_len}_wobble_gap_{gap}.pkl"
    shell:
        "python scripts/create_footprint_dict.py {input.footprint_capped} {output.footprint_capped_pkl}" 


rule build_footprint_matrix: 
    input:
        #tf_footprint_on_reads = "footprints/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_tf_footprint.bed",
        tf_footprint_on_reads = "footprint_min_size_selection/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_tf_footprint_min_length_10_wobble_gap_1.bed"
         
        #samflag_dict = "bismark_mapped/alignment_samflag_{sample}.pkl"
    params:
    output:
        methylation_matrix = "flank_footprint_matrix/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_methylation_matrix.tsv" ,
        methylation_matrix_real_footprint_len_pkl =\
              "flank_footprint_matrix/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_methylation_matrix_flen.pkl",
        methylation_matrix_footprint_vec_pkl = \
              "flank_footprint_matrix/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_methylation_matrix_footprint.pkl",
        strand_agnostic_footprint_vec_pkl = \
              "flank_footprint_matrix/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_methylation_matrix_footprint.strand.agnostic.pkl",
        strand_agnostic_footprint_matrix = \
              "flank_footprint_matrix/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_methylation_matrix.strand.agnostic.tsv", 
        footprint_length_at_bp_resolution_pkl = \
              "flank_footprint_matrix/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_methylation_matrix.footprint_len.bp.res.pkl", 
        footprint_length_at_bp_resolution_tsv = \
              "flank_footprint_matrix/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_methylation_matrix.footprint_len.bp.res.tsv", 
        #mvec_matrix = "flank_footprint_matrix/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_mvec_matrix.tsv"

    shell:
        #"python scripts/build_flank_methylation_matrix.py {input.tf_footprint_on_reads}"
        "python scripts/build_flank_methylation_matrix_v2.py {input.tf_footprint_on_reads}" 
        " {wildcards.lflank} {wildcards.rflank}"
        " {output.methylation_matrix}"
        " {output.methylation_matrix_real_footprint_len_pkl}"
        " {output.methylation_matrix_footprint_vec_pkl}"
        " {output.strand_agnostic_footprint_vec_pkl}"
        " {output.strand_agnostic_footprint_matrix}"
        " {output.footprint_length_at_bp_resolution_pkl}"
        " {output.footprint_length_at_bp_resolution_tsv}"
        #" {output.mvec_matrix}"

rule select_for_cluster:
    input:
        tf_footprint_on_reads = "footprints/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_tf_footprint.bed",
    params:
    output:
        #tf_footprint_in_cluster = "footprints/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_tf_footprint.bed",
        tf_footprint_in_cluster = "footprints/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_tf_footprint.cluster_{cl_id}.bed",
    shell:
       # "grep \"\^{wildcards.cl_id}\" {input.tf_footprint_on_reads} > {output.tf_footprint_in_cluster}" 
       "sh scripts/select_cluster.sh {input.tf_footprint_on_reads} {wildcards.cl_id} {output.tf_footprint_in_cluster}"
 
  
rule select_methylation_matrix_for_cluster: 
    input:
        methylation_matrix = "flank_footprint_matrix/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_methylation_matrix.tsv" 
        #tf_footprint_on_reads = "footprints/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_tf_footprint.bed",
    params:
    output:
        methylation_matrix_for_cl =  "flank_footprint_matrix/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_methylation_matrix.cluster_{cl_id}.tsv"
    shell:
        "sh scripts/select_cluster.sh {input.methylation_matrix} {wildcards.cl_id} {output.methylation_matrix_for_cl}"

    
rule cluster_footprint_vector_per_site:
    input: 
        #methylation_matrix = "flank_footprint_matrix/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix.tsv" 
        #methylation_matrix = "flank_footprint_matrix/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix.tsv" 
        methylation_matrix_for_cl =  "flank_footprint_matrix/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_methylation_matrix.cluster_{cl_id}.tsv"
    params:
        string_for_site = lambda wildcards: config["site_annotation"][wildcards.cl_id][wildcards.site] + "\^" + wildcards.cl_id
    output:
        site_specific_cluster_83_163 = "footprint_site_specific_clustering/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_83~163.tsv", 
        site_specific_cluster_99_147 = "footprint_site_specific_clustering/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_99~147.tsv"
    shell:
        "python scripts/clsuter_methylation_vector.py {input.methylation_matrix_for_cl}"
        " {params.string_for_site} {output.site_specific_cluster_83_163}"
        " {output.site_specific_cluster_99_147} {wildcards.nclust}"



def get_strand_info(wildcards): 
    peak_data = pd.read_csv(config["peak_annotation"], 
                            sep = "\s+", header = "infer")
    key = config["site_annotation"][wildcards.cl_id][wildcards.site] + "^" + wildcards.cl_id
    strand = peak_data.loc[peak_data.chr_loc == key,]["strand"].tolist()[0]
    return strand 
rule draw_footprint_vector_per_site:
    input:
        site_specific_cluster = "footprint_site_specific_clustering/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}.tsv",
    params:
        title =  lambda wildcards: "MNase Cluster " + wildcards.cl_id + " " +  config["site_annotation"][wildcards.cl_id][wildcards.site], 
        ylabel = lambda wildcards: config ["sam_flag_annotation"][wildcards.sam_flag], 
        strand = lambda wildcards: get_strand_info(wildcards)
    output:
        methylation_status_as_circles = "plots/footprint_methylation_status_for_clusters/methylation_status_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}.pdf", 
        average_methyaltion_all_in_one = "plots/footprint_average_methylation_status/avg_methylation_status_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}.pdf", 
        average_methyaltion_all_in_one_plot_tsv = "plots/footprint_average_methylation_status/avg_methylation_status_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}_plot.tsv", 
    shell:
        "Rscript scripts/draw_methylation_status_v2.R {input.site_specific_cluster}"
        " {output.methylation_status_as_circles} {wildcards.lflank} {wildcards.rflank}"
        " \"{params.title}\" \"{params.ylabel}\" {output.average_methyaltion_all_in_one}"
        " {output.average_methyaltion_all_in_one_plot_tsv} \"{params.strand}\""



rule plot_short_and_long_mnase_fragments_for_a_site: 
    input: 
        #short_mnase = "mnase_map/comb_50_to_reordered_20clust_147_50PNE_open_v1_peak_{ref_flank}_map.csv.gz_e.csv.gz", 
        #long_mnase = "mnase_map/comb_147_to_reordered_20clust_147_50PNE_open_v1_peak_{ref_flank}_map.csv.gz_e.csv.gz", 
        #short_mnase = "mnase_map/sm.sub.comb_50_c_to_reordered_20clust_147_50PNE_open_v1_peak_{ref_flank}_map.csv.gz_e.csv.gz", 
        #long_mnase = "mnase_map/sm.sub.comb_147_c_to_reordered_20clust_147_50PNE_open_v1_peak_{ref_flank}_map.csv.gz_e.csv.gz", 
        #short_mnase = "mnase_map/sm.sub.comb_50_c_to_all_peaks_{ref_flank}_map.csv.gz_e.csv.gz", 
        #long_mnase = "mnase_map/sm.sub.comb_147_c_to_all_peaks_{ref_flank}_map.csv.gz_e.csv.gz", 
        short_mnase = "mnase_map/sm.sub.comb_50_c_to_all_open_and_closed_mnase_peaks_{ref_flank}_map.csv.gz_e.csv.gz", 
        long_mnase = "mnase_map/sm.sub.comb_147_c_to_all_open_and_closed_mnase_peaks_{ref_flank}_map.csv.gz_e.csv.gz", 
        base_gnuplt_file = "utils/gnuplot_base_files/site_specific_mnase_eom.gplt"
    params:
        string_for_site = lambda wildcards: config["site_annotation"][wildcards.cl_id][wildcards.site] + "\^" + wildcards.cl_id, 
        title =  lambda wildcards: "MNase Cluster " + wildcards.cl_id + " " +  config["site_annotation"][wildcards.cl_id][wildcards.site]
      
    output:
        site_mnase_data = "site_specific_mnase_data/taken_from_{ref_flank}_cluster_{cl_id}_span_from_center_{span}_vline_left_{left}_right_{right}_site_{site}.tsv",
        both_in_one_pdf =  "plots/site_specific_mnase/taken_from_{ref_flank}_cluster_{cl_id}_span_from_center_{span}_vline_left_{left}_right_{right}_site_{site}.pdf",
        both_in_one_gplt = "plots/site_specific_mnase/taken_from_{ref_flank}_cluster_{cl_id}_span_from_center_{span}_vline_left_{left}_right_{right}_site_{site}.gplt",
        mnase_top_panel = "plots/site_specific_mnase/taken_from_{ref_flank}_cluster_{cl_id}_span_from_center_{span}_vline_left_{left}_right_{right}_site_{site}.top.panel.pdf"
    shell:
        "sh scripts/plot_site_specific_mnase.sh {input.short_mnase} {input.long_mnase} {input.base_gnuplt_file}"
        " {output.site_mnase_data} {output.both_in_one_gplt} {output.both_in_one_pdf} \"{params.string_for_site}\""
        " {wildcards.span} {wildcards.left} {wildcards.right} \"'{params.title}'\" {output.mnase_top_panel}"


rule plot_methylation_and_mnase_together:
    input: 
        #average_methyaltion_all_in_one_plot_tsv = "plots/average_methylation_status/avg_methylation_status_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}_plot.tsv", 
        average_methyaltion_all_in_one_plot_tsv = "plots/footprint_average_methylation_status/avg_methylation_status_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}_plot.tsv", 
        site_mnase_data = "site_specific_mnase_data/taken_from_{ref_flank}_cluster_{cl_id}_span_from_center_{span}_vline_left_{left}_right_{right}_site_{site}.tsv",
        
    params:
    output:
        mnase_and_methylation_pdf = "plots/mnase_and_avg_methylation/methylation_{sample}_to_{bed}_cl_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_nclust_{nclust}_sflag_{sam_flag}_mnase_{ref_flank}_span_{span}_mlf_{left}_mrf_{right}_grey_box_left_{left_x}_right_{right_x}.pdf", # left_x: comma separated x values for left x-coordinate; similar for right_x ; mlf: mnase left vline ; 
        mnase_and_methylation_png = "plots/mnase_and_avg_methylation/methylation_{sample}_to_{bed}_cl_{cl_id}_lf_{lflank}_rf_{rflank}_site_{site}_nclust_{nclust}_sflag_{sam_flag}_mnase_{ref_flank}_span_{span}_mlf_{left}_mrf_{right}_grey_box_left_{left_x}_right_{right_x}.png" 
    shell:
        "Rscript scripts/plot_mnase_and_methylation.R" 
        " {input.average_methyaltion_all_in_one_plot_tsv}"
        " {input.site_mnase_data} {output.mnase_and_methylation_pdf}"
        " {output.mnase_and_methylation_png}"
        " \"{wildcards.left_x}\" \"{wildcards.right_x}\""

rule draw_footprint_vector_per_site_with_mnase:
    input:
        site_specific_cluster = "footprint_site_specific_clustering/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}.tsv",
        site_mnase_data = "site_specific_mnase_data/taken_from_{ref_flank}_cluster_{cl_id}_span_from_center_{span}_vline_left_{left}_right_{right}_site_{site}.tsv",
        
    params:
        title =  lambda wildcards: "MNase Cluster " + wildcards.cl_id + " " +  config["site_annotation"][wildcards.cl_id][wildcards.site], 
        ylabel = lambda wildcards: config ["sam_flag_annotation"][wildcards.sam_flag], 
        strand = lambda wildcards: get_strand_info(wildcards)
    output:
        methylation_status_as_circles = "plots/mnase_and_footprint_methylation_status_for_clusters/methylation_status_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}_mnase_ref_{ref_flank}_span_{span}_vline_l_{left}_r_{right}.pdf", 
        average_methyaltion_all_in_one = "plots/mnase_and_footprint_average_methylation_status/avg_methylation_status_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}_mnase_ref_{ref_flank}_span_{span}_vline_l_{left}_r_{right}.pdf", 
        average_methyaltion_all_in_one_plot_tsv = "plots/footprint_average_methylation_status/avg_methylation_status_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}_mnase_ref_{ref_flank}_span_{span}_vline_l_{left}_r_{right}_plot.tsv", 
    shell:
        "Rscript scripts/draw_methylation_status_v3.R {input.site_specific_cluster}"
        " {output.methylation_status_as_circles} {wildcards.lflank} {wildcards.rflank}"
        " \"{params.title}\" \"{params.ylabel}\" {output.average_methyaltion_all_in_one}"
        " {output.average_methyaltion_all_in_one_plot_tsv} \"{params.strand}\""
        " {input.site_mnase_data} {wildcards.left} {wildcards.right} {wildcards.site}"

rule preapre_extended_footprint_tsv: 
    input:
        site_specific_cluster = "footprint_site_specific_clustering/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}.tsv", 
        #extended_matrix_pkl = "flank_footprint_matrix/{sample}_to_{bed}_lf_{lextend}_rf_{rextend}_methylation_matrix_footprint.pkl"
        extended_matrix_pkl = "fragments_mapped_to_enhancer_center/{sample}_to_{bed}_footprint.pkl"
    params:
        strand = lambda wildcards: get_strand_info(wildcards)
    output:
        extended_footprint_tsv = "footprint_site_specific_clustering/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}_extend_from_peak_left_{lextend}_right_{rextend}.tsv" ,
        extended_footprint_tsv_verbose = "footprint_site_specific_clustering/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}_extend_from_peak_left_{lextend}_right_{rextend}.verbose.tsv" 
    shell:
        "python scripts/extend_footprint.py {input.site_specific_cluster}"
        " {input.extended_matrix_pkl} {wildcards.lflank} {wildcards.rflank}"
        " {wildcards.lextend} {wildcards.rextend} {output.extended_footprint_tsv}"
        " {params.strand} {output.extended_footprint_tsv_verbose}" 
         
rule draw_footprint_vector_per_site_with_mnase_extended:
    input:
        extended_footprint_tsv = "footprint_site_specific_clustering/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}_extend_from_peak_left_{lextend}_right_{rextend}.tsv", 
        site_mnase_data = "site_specific_mnase_data/taken_from_{ref_flank}_cluster_{cl_id}_span_from_center_{span}_vline_left_{left}_right_{right}_site_{site}.tsv",
        
    params:
        title =  lambda wildcards: "MNase Cluster " + wildcards.cl_id + " " +  config["site_annotation"][wildcards.cl_id][wildcards.site], 
        ylabel = lambda wildcards: config ["sam_flag_annotation"][wildcards.sam_flag], 
        strand = lambda wildcards: get_strand_info(wildcards)
    output:
        methylation_status_as_circles = "plots/extended_mnase_and_footprint_methylation_status_for_clusters/methylation_status_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}_extend_from_peak_left_{lextend}_right_{rextend}_mnase_ref_{ref_flank}_span_{span}_vline_l_{left}_r_{right}.extended.pdf", 
        average_methyaltion_all_in_one = "plots/extended_mnase_and_footprint_average_methylation_status/avg_methylation_status_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}_extend_from_peak_left_{lextend}_right_{rextend}_mnase_ref_{ref_flank}_span_{span}_vline_l_{left}_r_{right}.extended.pdf", 
        average_methyaltion_all_in_one_plot_tsv = "plots/extended_footprint_average_methylation_status/avg_methylation_status_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}_extend_from_peak_left_{lextend}_right_{rextend}_mnase_ref_{ref_flank}_span_{span}_vline_l_{left}_r_{right}_plot.extended.tsv",
        
    shell:
        "Rscript scripts/draw_methylation_status_v3_extend.R {input.extended_footprint_tsv}"
#        " {output.methylation_status_as_circles} {wildcards.lflank} {wildcards.rflank}"
        " {output.methylation_status_as_circles} {wildcards.lextend} {wildcards.rextend}"
        " \"{params.title}\" \"{params.ylabel}\" {output.average_methyaltion_all_in_one}"
        " {output.average_methyaltion_all_in_one_plot_tsv} \"{params.strand}\""
        " {input.site_mnase_data} {wildcards.lflank} {wildcards.rflank}"
        " {wildcards.site} {wildcards.sam_flag}"
          

rule draw_footprint_vector_per_site_with_mnase_with_extended_flank:
    input:
        site_specific_cluster = "footprint_site_specific_clustering/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}.tsv",
        site_mnase_data = "site_specific_mnase_data/taken_from_{ref_flank}_cluster_{cl_id}_span_from_center_{span}_vline_left_{left}_right_{right}_site_{site}.tsv",
        
    params:
        title =  lambda wildcards: "MNase Cluster " + wildcards.cl_id + " " +  config["site_annotation"][wildcards.cl_id][wildcards.site], 
        ylabel = lambda wildcards: config ["sam_flag_annotation"][wildcards.sam_flag], 
        strand = lambda wildcards: get_strand_info(wildcards)
    output:
        methylation_status_as_circles = "plots/mnase_and_footprint_methylation_status_for_clusters/methylation_status_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}_mnase_ref_{ref_flank}_span_{span}_vline_l_{left}_r_{right}.pdf", 
        average_methyaltion_all_in_one = "plots/mnase_and_footprint_average_methylation_status/avg_methylation_status_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}_mnase_ref_{ref_flank}_span_{span}_vline_l_{left}_r_{right}.pdf", 
        average_methyaltion_all_in_one_plot_tsv = "plots/footprint_average_methylation_status/avg_methylation_status_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}_mnase_ref_{ref_flank}_span_{span}_vline_l_{left}_r_{right}_plot.tsv", 
    shell:
        "Rscript scripts/draw_methylation_status_v3.R {input.site_specific_cluster}"
        " {output.methylation_status_as_circles} {wildcards.lflank} {wildcards.rflank}"
        " \"{params.title}\" \"{params.ylabel}\" {output.average_methyaltion_all_in_one}"
        " {output.average_methyaltion_all_in_one_plot_tsv} \"{params.strand}\""
        " {input.site_mnase_data} {wildcards.left} {wildcards.right} {wildcards.site}"




# Get the true length profile with kde

rule real_footprint_lengths:
     input:
         site_specific_cluster = "footprint_site_specific_clustering/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}.tsv",
         #footprint_count_file = "footprints/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_tf_footprint.count" : for scripts/get_real_footprint_lengths_from_clustered_reads.py
         methylation_matrix_real_footprint_len_pkl =\
              "flank_footprint_matrix/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_methylation_matrix_flen.pkl" 
     params:
     output: 
         real_footprint_length_tsv = "actual_footprint_length_in_clusters/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}_footprint_length.tsv", 
         real_footprint_length_and_per_orange = "actual_footprint_length_in_clusters/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}_footprint_length.per_orange.tsv"
     shell: 
         #"python scripts/get_real_footprint_lengths_from_clustered_reads.py"
         #" {input.site_specific_cluster} {input.footprint_count_file}"
         #" \"-{wildcards.lflank}\" \"{wildcards.rflank}\""
         #" {output.real_footprint_length_tsv}"
         "python scripts/get_real_footprint_lengths_from_clustered_reads_v2.py"
         " {input.site_specific_cluster} {input.methylation_matrix_real_footprint_len_pkl}"
         " {output.real_footprint_length_tsv}"
         " {output.real_footprint_length_and_per_orange}"
rule plot_footprint_length_vs_per_orange:
    input:
         real_footprint_length_and_per_orange = "actual_footprint_length_in_clusters/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}_footprint_length.per_orange.tsv"
    params:
    output:
        length_vs_per_orange_pdf =  "plots/footprint_length_vs_per_orange/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}_footprint_length.per_orange.pdf",
        length_vs_per_orange_png =  "plots/footprint_length_vs_per_orange/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}_footprint_length.per_orange.png",
    shell:
        "Rscript scripts/footprint_length_vs_per_orange.R"
        " {input.real_footprint_length_and_per_orange}"
        " {output.length_vs_per_orange_pdf} {output.length_vs_per_orange_png}"
    

rule kde_real_footprint_length_distribution:
    input:
        real_footprint_length_tsv = "actual_footprint_length_in_clusters/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}_footprint_length.tsv", 
        col_data = "utils/gnuplot_base_files/gnuplot_color_selection.list", 
        line_type_list = "utils/r_base_files/line_style_list.list"
    params:
        bw = 10,
        n = 40, 
        start = 0, 
        end = 200
    output:
        kde_dist_individual_pdf = "plots/actual_footprint_length_in_clusters/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}_footprint_length_kde.ind.pdf",
        kde_dist_individual_png = "plots/actual_footprint_length_in_clusters/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}_footprint_length_kde.ind.png",  
        kde_dist_combined_pdf = "plots/actual_footprint_length_in_clusters/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}_footprint_length_kde.comb.pdf",
        kde_dist_combined_png = "plots/actual_footprint_length_in_clusters/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}_footprint_length_kde.comb.png", 
        kde_df_tsv = "plots/actual_footprint_length_in_clusters/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}_footprint_length_kde.comb.tsv" , 
        hist_of_all_lengths_combined_pdf = "plots/actual_footprint_length_in_clusters/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}_footprint_length_hist.all.pdf",
        hist_of_all_lengths_combined_png = "plots/actual_footprint_length_in_clusters/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}_footprint_length_hist.all.png", 
        hist_of_all_lengths_cap_200_combined_pdf = "plots/actual_footprint_length_in_clusters/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}_footprint_length_hist.cap_200.pdf",
        hist_of_all_lengths_cap_200_combined_png = "plots/actual_footprint_length_in_clusters/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}_footprint_length_hist.cap_200.png", 
        
        
    shell:
        "Rscript scripts/kde_real_footprint_length.R {input.real_footprint_length_tsv}"
        " \"1\" \"2\" {params.bw} {input.col_data} \"Cluster\""
        " \"Footprint Length[bp]\" \"Density [A.U.]\""
        " {input.line_type_list} {output.kde_dist_individual_png}"
        " {output.kde_dist_individual_pdf}"
        " {wildcards.nclust} {output.kde_dist_combined_png}"
        " {output.kde_dist_combined_pdf} {output.kde_df_tsv}"
        " {output.hist_of_all_lengths_combined_pdf} {output.hist_of_all_lengths_combined_png}" 
        " {wildcards.site}"
        " {output.hist_of_all_lengths_cap_200_combined_pdf} {output.hist_of_all_lengths_cap_200_combined_png}"

# python scripts/calculate_percentage_tf_nuc_and_naked_dna.py footprint_site_specific_clustering/merged_S2_to_reordered_20clust_147_50PNE_open_v1_peak_cluster_12_lf_100_rf_50_methylation_matrix_for_site_peak_2299_nclust_3_sam_flag_99~147.tsv plots/actual_footprint_length_in_clusters/merged_S2_to_reordered_20clust_147_50PNE_open_v1_peak_cluster_12_lf_100_rf_50_methylation_matrix_for_site_peak_2299_nclust_3_sam_flag_99~147_footprint_length_kde.comb.tsv tmp/wq.json tmp/wq.tsv tmp/per.tsv
rule calculate_percentage_footprints: 
    input:
        site_specific_cluster = "footprint_site_specific_clustering/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}.tsv",
        kde_df_tsv = "plots/actual_footprint_length_in_clusters/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}_footprint_length_kde.comb.tsv"
    params:
    output:
        data_json_file = "auto_cluster_assingment/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}.for_auto.json",
        data_tsv_file = "auto_cluster_assingment/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}.for_auto.tsv",
        percentage_tsv = "auto_cluster_assingment/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}.auto_per.tsv"
    
    shell:
        "python scripts/calculate_percentage_tf_nuc_and_naked_dna.py"
        " {input.site_specific_cluster} {input.kde_df_tsv}"
        " {output.data_json_file} {output.data_tsv_file} {output.percentage_tsv}"
        " {wildcards.site}"

rule calculate_percentage_orange_per_read_in_a_peak:
    input:
        site_specific_cluster = "footprint_site_specific_clustering/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}.tsv",
    params:
    output:
        per_orange_file = "percetage_orange_per_read/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}.per_orange.tsv"
    shell:
        "python scripts/percentage_orange_per_read.py {input.site_specific_cluster} {output.per_orange_file}"

rule plot_percentage_orange_for_a_peak: 
    input:
        per_orange_file = "percetage_orange_per_read/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}.per_orange.tsv"
    params: 
    output:
        per_orange_hist_pdf = "plots/percentage_orange_plots/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}.per_orange.hist.pdf", 
        per_orange_hist_png = "plots/percentage_orange_plots/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}.per_orange.hist.png", 
    shell:
        "Rscript scripts/hist_percentage_orange.R {input.per_orange_file}"
        " {output.per_orange_hist_pdf} {output.per_orange_hist_png}"
rule merge_different_plots: 
    input:
        per_orange_hist_png_99= "plots/percentage_orange_plots/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_99~147.per_orange.hist.png", 
        per_orange_hist_png_83 = "plots/percentage_orange_plots/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_83~163.per_orange.hist.png", 
        hist_of_all_lengths_cap_200_combined_png_99 = "plots/actual_footprint_length_in_clusters/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_99~147_footprint_length_hist.cap_200.png", 
        hist_of_all_lengths_cap_200_combined_png_83 = "plots/actual_footprint_length_in_clusters/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_83~163_footprint_length_hist.cap_200.png", 
        dot_plot_99  = "plots/extended_mnase_and_footprint_methylation_status_for_clusters/methylation_status_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_99~147_extend_from_peak_left_{lextend}_right_{rextend}_mnase_ref_{ref_flank}_span_{span}_vline_l_{left}_r_{right}.extended.pdf", 
        dot_plot_83 = "plots/extended_mnase_and_footprint_methylation_status_for_clusters/methylation_status_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_83~163_extend_from_peak_left_{lextend}_right_{rextend}_mnase_ref_{ref_flank}_span_{span}_vline_l_{left}_r_{right}.extended.pdf", 
        length_vs_per_orange_png_99 =  "plots/footprint_length_vs_per_orange/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_99~147_footprint_length.per_orange.png",
        length_vs_per_orange_png_83 =  "plots/footprint_length_vs_per_orange/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_83~163_footprint_length.per_orange.png",
        
    params:
    output:
        combined_plots = "plots/combined_plots/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_extend_from_peak_left_{lextend}_right_{rextend}_mnase_ref_{ref_flank}_span_{span}_vline_l_{left}_r_{right}.combined.png" 
    shell:
        "sh scripts/combine_plots.sh {input.per_orange_hist_png_99} {input.per_orange_hist_png_83} "
        " {input.hist_of_all_lengths_cap_200_combined_png_99} {input.hist_of_all_lengths_cap_200_combined_png_83}"  
        " {input.dot_plot_99} {input.dot_plot_83}"
        " {output.combined_plots}"
        " {input.length_vs_per_orange_png_99} {input.length_vs_per_orange_png_83}"
  
############# implement inverted U shpape algorithm ###################
rule footprint_len_per_bp_for_a_site: 
    input:
        site_specific_cluster = "footprint_site_specific_clustering/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}.tsv",
        footprint_length_at_bp_resolution_pkl =  "flank_footprint_matrix/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_methylation_matrix.footprint_len.bp.res.pkl",
    params:
    output:
        footprint_len_per_bp =  "footprint_length_per_bp_for_sites/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}.footprint_per_bp.tsv"
    shell:
        "python scripts/prepare_data_for_inverted_u.py {input.site_specific_cluster}"
        " {input.footprint_length_at_bp_resolution_pkl} {output.footprint_len_per_bp}"
rule plot_inverse_u:
    input: 
        footprint_len_per_bp =  "footprint_length_per_bp_for_sites/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}.footprint_per_bp.tsv",
    params:
    output:
        mean_and_median_u_pdf = "plots/inverse_u_plots/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}.inverse_u.mean_and_median.pdf",
        mean_and_median_u_png = "plots/inverse_u_plots/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}.inverse_u.mean_and_median.png",
    shell:
        "Rscript scripts/plot_inverse_u.R {input.footprint_len_per_bp}"
        " {output.mean_and_median_u_pdf}"
        " {wildcards.lflank} {wildcards.rflank} {wildcards.nclust}"
        " {output.mean_and_median_u_png}"        


rule combine_inverse_u_and_dot_plot: 
    input:
        mean_and_median_u_png_99 = "plots/inverse_u_plots/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_99~147.inverse_u.mean_and_median.png",
        mean_and_median_u_png_83 = "plots/inverse_u_plots/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_83~163.inverse_u.mean_and_median.png",
        dot_plot_99  = "plots/extended_mnase_and_footprint_methylation_status_for_clusters/methylation_status_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_99~147_extend_from_peak_left_{lextend}_right_{rextend}_mnase_ref_{ref_flank}_span_{span}_vline_l_{left}_r_{right}.extended.pdf", 
        dot_plot_83 = "plots/extended_mnase_and_footprint_methylation_status_for_clusters/methylation_status_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_83~163_extend_from_peak_left_{lextend}_right_{rextend}_mnase_ref_{ref_flank}_span_{span}_vline_l_{left}_r_{right}.extended.pdf", 
    
    params:
    output:
        combined_plot = "plots/combined_u_and_dot_plots/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_extend_from_peak_left_{lextend}_right_{rextend}_mnase_ref_{ref_flank}_span_{span}_vline_l_{left}_r_{right}.dot_and_u.png"
    shell:
        "sh scripts/combine_u_and_dot_plots.sh {input.mean_and_median_u_png_99}"
        " {input.mean_and_median_u_png_83} {input.dot_plot_99} {input.dot_plot_83}"
        " {output.combined_plot}"

        

##########################################################
#
# Aoid below
#
########################################################## 

#rule build_flank_matrix: 
#    input:
#        selected_flank = "flank_with_methyation_status/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_mapped.bed", 
#        samflag_dict = "bismark_mapped/alignment_samflag_{sample}.pkl"
#    params:
#    output:
#        methylation_matrix = "flank_methylation_matrix/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix.tsv"
#    shell:
#        "python scripts/build_flank_methylation_matrix.py {input.selected_flank}"
#        " {input.samflag_dict} {wildcards.lflank} {wildcards.rflank}"
#        " {output.methylation_matrix}"
#
#
#rule clsuter_methylation_vector_per_site:
#    input:
#        methylation_matrix = "flank_methylation_matrix/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix.tsv"
#    params: 
#        string_for_site = lambda wildcards: config["site_annotation"][wildcards.cl_id][wildcards.site] + "\^" + wildcards.cl_id
#    output: 
#        site_specific_cluster_83_163 = "site_specific_clustering/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_83~163.tsv", 
#        site_specific_cluster_99_147 = "site_specific_clustering/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_99~147.tsv"
#
#    shell:
#        "python scripts/clsuter_methylation_vector.py {input.methylation_matrix}"
#        " {params.string_for_site} {output.site_specific_cluster_83_163} "
#        " {output.site_specific_cluster_99_147} {wildcards.nclust}"
#rule draw_read_methylation_status:
#    input:
#        site_specific_cluster = "site_specific_clustering/{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}.tsv",
#    params:
#        title =  lambda wildcards: "MNase Cluster " + wildcards.cl_id + " " +  config["site_annotation"][wildcards.cl_id][wildcards.site], 
#        ylabel = lambda wildcards: config ["sam_flag_annotation"][wildcards.sam_flag]
#    output:
#        methylation_status_as_circles = "plots/methylation_status_for_clusters/methylation_status_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}.pdf", 
#        average_methyaltion_all_in_one = "plots/average_methylation_status/avg_methylation_status_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}.pdf", 
#        average_methyaltion_all_in_one_plot_tsv = "plots/average_methylation_status/avg_methylation_status_{sample}_to_{bed}_cluster_{cl_id}_lf_{lflank}_rf_{rflank}_methylation_matrix_for_site_{site}_nclust_{nclust}_sam_flag_{sam_flag}_plot.tsv", 
#    shell:
#        "Rscript scripts/draw_methylation_status_v2.R {input.site_specific_cluster}"
#        " {output.methylation_status_as_circles} {wildcards.lflank} {wildcards.rflank}"
#        " \"{params.title}\" \"{params.ylabel}\" {output.average_methyaltion_all_in_one}"
#        " {output.average_methyaltion_all_in_one_plot_tsv}"
#
#        
