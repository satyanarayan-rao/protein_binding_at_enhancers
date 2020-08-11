rule get_enhancer_mapped_reads_bed:
    input:
        tf_footprint_on_reads = "footprint_min_size_selection/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_tf_footprint_min_length_10_wobble_gap_1.bed"
    params:
    output:
        reads_bed = "obs_vs_theoretical_analysis/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_reads.bed"
    shell:
        "sh scripts/prepare_mapped_reads_bed.sh {input.tf_footprint_on_reads}"
        " {output.reads_bed}"

rule get_read_genomic_sequence_in_tsv:
    input:
        reads_bed = "obs_vs_theoretical_analysis/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_reads.bed",
        dm3_fasta = "/beevol/home/satyanarr/workplace/data/ucsc/dm/dm3/dm3.fa"
    params:
    output:
        genomic_seq_tsv = "obs_vs_theoretical_analysis/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_read_genomic_seq.tsv"
    shell:
        "bedtools getfasta -fi {input.dm3_fasta} -fo - -bed {input.reads_bed} -tab -name+ | sed 's/::/\t/g;s/:/\t/g;s/-/\t/g' | awk '{{print $2\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$1}}' > {output.genomic_seq_tsv}"        

rule get_theoretical_footprints:
    input:
        genomic_seq_tsv = "obs_vs_theoretical_analysis/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_read_genomic_seq.tsv",
    params:
    output:
        theoretical_footprints = "obs_vs_theoretical_analysis/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_theoretical_footprints.tsv"
    shell:
        "python scripts/footprint_len_in_genomic_seq_bed_coordinates.py"
        " {input.genomic_seq_tsv} {output.theoretical_footprints}" 

rule prepare_enhancer_peak_flank_bed:
    input:
        input_bed = "input_bed/{bed}.bed", 
        dm3_chrom_size = "metadata/dm3.chrom.sizes"
    params:
    output:
        flanked_bed = "input_bed/{bed}_lf_{lflank}_rf_{rflank}.bed" 
    shell:
        "bedtools slop -l {wildcards.lflank} -r {wildcards.rflank} -s -g {input.dm3_chrom_size} -i {input.input_bed} > {output.flanked_bed}"   

rule intesect_theoretical_footprint_to_peak:
    input: 
        theoretical_footprints = "obs_vs_theoretical_analysis/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_theoretical_footprints.tsv", 
        flanked_bed = "input_bed/{bed}_lf_{lflank}_rf_{rflank}.bed" 
    params:
    output:
        intersected_out = "obs_vs_theoretical_analysis/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_read_intersect.bed" 
    shell:
        "bedtools intersect -a {input.theoretical_footprints} -b {input.flanked_bed} -wa -wb > {output.intersected_out}" 

rule observed_footprints_on_read:
    input:
        tf_footprint_on_reads = "footprint_min_size_selection/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_tf_footprint_min_length_10_wobble_gap_1.bed"
    params:
    output:
        obs_footprints = "obs_vs_theoretical_analysis/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_observed_footprints.bed"
    shell:
        #"python scripts/get_length_and_loc_footprint.py {input.tf_footprint_on_reads} {output.obs_footprints}"
        "python scripts/get_length_and_loc_footprint_include_occluded.py {input.tf_footprint_on_reads} {output.obs_footprints}"

rule intersect_oserved_footprints_to_peak:
    input:
        obs_footprints = "obs_vs_theoretical_analysis/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_observed_footprints.bed",
        flanked_bed = "input_bed/{bed}_lf_{lflank}_rf_{rflank}.bed" 
        
    params:
    output:
        intersected_out = "obs_vs_theoretical_analysis/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_footprint_intersect.bed"
    shell: 
        "bedtools intersect -a {input.obs_footprints} -b {input.flanked_bed} -wa -wb > {output.intersected_out}" 

rule prepare_observed_vs_expected_file:
    input:
        theoretical = "obs_vs_theoretical_analysis/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_read_intersect.bed",
        observed = "obs_vs_theoretical_analysis/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_footprint_intersect.bed"
    params:
    output:
        obs_vs_exp_file = "obs_vs_theoretical_analysis/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_obs_vs_expected.tsv"
    shell:
        "sh scripts/obs_vs_theoretical_footprints.sh {input.theoretical} {input.observed}"
        " {output.obs_vs_exp_file}"

rule count_footprint_length:
    input:
        intersect_file = "obs_vs_theoretical_analysis/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_{ftype}_intersect.bed"
    params:
        column_id = 4, 
        max_lim = 350 
    output:
        count_file = "obs_vs_theoretical_analysis/{sample}_to_{bed}_lf_{lflank}_rf_{rflank}_{ftype}_intersect.cnt" 
    shell:
        "python scripts/frequency_of_column.py {input.intersect_file}"
        " {output.count_file} {params.column_id} {params.max_lim} {wildcards.ftype} {wildcards.bed}"
        

########################### Do it on the whole enhancers from open and closed ####################################
# Motivating behind doing this separately is that with current implementation I am limited to just peak center   #
#                                                                                                                #
# And also I have to process data in compressed form as theoretical footprint files will be very large           #
##################################################################################################################


rule get_whole_enh_mapped_reads_in_bed:
    input:
        footprint_bed = "fragments_mapped_to_enhancer_center/{sample}_to_{bed}_footprint.bed",  
    params:
    output:
        reads_bed = "obs_vs_theoretical_analysis/{sample}_to_{bed}_reads.bed"
    shell:
        "sh scripts/prepare_mapped_reads_bed.sh {input.footprint_bed}"
        " {output.reads_bed}"
        
rule get_read_genomic_sequence_in_tsv_whole_enh:
    input:
        reads_bed = "obs_vs_theoretical_analysis/{sample}_to_{bed}_reads.bed",
        dm3_fasta = "/beevol/home/satyanarr/workplace/data/ucsc/dm/dm3/dm3.fa"
    params:
    output:
        genomic_seq_tsv = "obs_vs_theoretical_analysis/{sample}_to_{bed}_read_genomic_seq.tsv"
    shell:
        "bedtools getfasta -fi {input.dm3_fasta} -fo - -bed {input.reads_bed} -tab -name+ | sed 's/::/\t/g;s/:/\t/g;s/-/\t/g' | awk '{{print $2\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$1}}' > {output.genomic_seq_tsv}"        



rule get_theoretical_footprints_whole_enh:
    input:
        genomic_seq_tsv = "obs_vs_theoretical_analysis/{sample}_to_{bed}_read_genomic_seq.tsv",
    params:
    output:
        theoretical_footprints = "obs_vs_theoretical_analysis/{sample}_to_{bed}_theoretical_footprints.tsv.gz"
    shell:
        "python scripts/footprint_len_in_genomic_seq_bed_coordinates_gzip.py"
        " {input.genomic_seq_tsv} {output.theoretical_footprints}" 

rule freq_of_theoretical_footprint_whole_enh:
    input:
        length_file = "obs_vs_theoretical_analysis/{sample}_to_{bed}_theoretical_footprints.tsv.gz" 
    params:
        column_id = 4, 
        max_lim = 350 
    output:
        freq_file = "obs_vs_theoretical_analysis/{sample}_to_{bed}_theoretical_footprints.freq"
    shell:
        "python scripts/frequency_of_column_gzip.py {input.length_file}"
        " {output.freq_file} {params.column_id} {params.max_lim} \"theoretical\" {wildcards.bed}"



############################### Theoretical ends here ############################


############################### Observed start here ##############################


rule fix_wobble_whole_enh:
    input:
        footprint_bed = "fragments_mapped_to_enhancer_center/{sample}_to_{bed}_footprint.bed"
    params:
    output:
        wobble_fixed = "fragments_mapped_to_enhancer_center/{sample}_to_{bed}_wobble_fixed_footprint_gap_{gap}.bed"
    shell:
        "python scripts/fix_wobble.py {input.footprint_bed}"
        " {wildcards.gap} {output.wobble_fixed}" 

rule min_cap_footprint_whole_enh:
    input:
        wobble_fixed = "fragments_mapped_to_enhancer_center/{sample}_to_{bed}_wobble_fixed_footprint_gap_{gap}.bed"
    output:
        footprint_capped = "fragments_mapped_to_enhancer_center/{sample}_to_{bed}_footprint_min_length_{min_len}_wobble_gap_{gap}.bed"    
    shell:
        "python scripts/keep_all_footprint_ge_th.py {input.wobble_fixed}"
        " {wildcards.min_len} {output.footprint_capped}"


rule observed_footprints_on_read_whole_enh:
    input:
        all_footprints = "fragments_mapped_to_enhancer_center/{sample}_to_{bed}_footprint_min_length_{min_len}_wobble_gap_{gap}.bed"
    params:
    output:
        length_file = "fragments_mapped_to_enhancer_center/{sample}_to_{bed}_footprint_min_length_{min_len}_wobble_gap_{gap}.cnt" 
    shell:
        "python scripts/get_length_and_loc_footprint_include_occluded.py {input.all_footprints} {output.length_file}"

rule freq_of_footprint_whole_enh:
    input:
        length_file = "fragments_mapped_to_enhancer_center/{sample}_to_{bed}_footprint_min_length_{min_len}_wobble_gap_{gap}.cnt" 
    params:
        column_id = 4, 
        max_lim = 350 
    output:
        freq_file = "fragments_mapped_to_enhancer_center/{sample}_to_{bed}_footprint_min_length_{min_len}_wobble_gap_{gap}.freq"
    shell:
        "python scripts/frequency_of_column.py {input.length_file}"
        " {output.freq_file} {params.column_id} {params.max_lim} \"footprint\" {wildcards.bed}"


#rule plot_obs_vs_theoretical_whole_enh:
#    input:
#        closed_footprint_freq_file = "fragments_mapped_to_enhancer_center/suppressed_merged_S2_to_closed_S2_STARRSEQ_dummy_cl_footprint_min_length_10_wobble_gap_1.freq", 
#        open_footprint_freq_file = "fragments_mapped_to_enhancer_center/suppressed_merged_S2_to_open_S2_STARRSEQ_dummy_cl_footprint_min_length_10_wobble_gap_1.freq", 
#        closed_theoretical_freq_file = "obs_vs_theoretical_analysis/suppressed_merged_S2_to_closed_S2_STARRSEQ_dummy_cl_theoretical_footprints.freq", 
#        open_theoretical_freq_file = "obs_vs_theoretical_analysis/suppressed_merged_S2_to_open_S2_STARRSEQ_dummy_cl_theoretical_footprints.freq", 
#        
#    params:
#        labels = "Obs-Closed@Obs-Open@Theoretical-Closed@Theoretical-Open"
#    output:
#        obs_vs_theoretical_whole_enh_pdf = "plots/obs_vs_theoretical/whole_enh.pdf",
#        obs_vs_theoretical_whole_enh_png = "plots/obs_vs_theoretical/whole_enh.png",
#        obs_vs_theoretical_whole_enh_tsv = "plots/obs_vs_theoretical/whole_enh.tsv",
#    shell:
#        "Rscript scripts/plot_obs_vs_theoretical.R "
#        " \"{input.closed_footprint_freq_file} {input.open_footprint_freq_file} {input.closed_theoretical_freq_file} {input.open_theoretical_freq_file}\""  
#        " {params.labels} {output.obs_vs_theoretical_whole_enh_pdf}"
#        " {output.obs_vs_theoretical_whole_enh_png}"
#        " {output.obs_vs_theoretical_whole_enh_tsv}"  
    


#rule prepare_obs_vs_theory_tsv_whole_enh:
#    input:
#        obs_closed = "fragments_mapped_to_enhancer_center/suppressed_merged_S2_to_closed_S2_STARRSEQ_dummy_cl_footprint_min_length_10_wobble_gap_1.cnt",
#        obs_open = "fragments_mapped_to_enhancer_center/suppressed_merged_S2_to_open_S2_STARRSEQ_dummy_cl_footprint_min_length_10_wobble_gap_1.cnt",
#        theoretical_closed = "obs_vs_theoretical_analysis/suppressed_merged_S2_to_closed_S2_STARRSEQ_dummy_cl_theoretical_footprints.tsv.gz", 
#        theoretical_open = "obs_vs_theoretical_analysis/suppressed_merged_S2_to_open_S2_STARRSEQ_dummy_cl_theoretical_footprints.tsv.gz",  
#    params:
#        labels = "Obs-Closed@obs_open@Theoretical-Closed@Theoretical-Open"
#    output:
#        obs_vs_theo_whole_enh = "obs_vs_theoretical_analysis/obs_vs_theoretical_whole_enh.tsv"
#    shell:
#        "python scripts/prepare_obs_vs_theory_tsv_whole_enh.py"
#        " \"{input.obs_closed} {input.obs_open} {input.theoretical_closed} {input.theoretical_open}\"" 
#        "{params.labels}"
        

#rule count_footprint_length_whole_enh:
#    input:
#    params:
#        column_id = 4, 
#        max_lim = 350
#    output:
#        count_file = "fragments_mapped_to_enhancer_center/{sample}_to_{bed}_footprint_min_length_{min_len}_wobble_gap_{gap}.cnt" 
#    shell:
#        "python scripts/frequency_of_column.py {input.all_footprints}"
#        " {output.count_file} {params.column_id} {params.max_lim} \"footprint\" {wildcards.bed}"
#
