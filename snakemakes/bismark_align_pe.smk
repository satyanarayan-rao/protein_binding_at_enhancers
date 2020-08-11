rule bismark_align_pe:
    input:
        read1 = "trimmed/{sample}_R1_val_1.fq.gz", 
        read2 = "trimmed/{sample}_R2_val_2.fq.gz"
    params:
        genome = "/beevol/home/satyanarr/workplace/data/ucsc/dm/dm3/", 
        out_dir = "bismark_mapped",
        extra = "--gzip --no_dovetail -p 4",
        basename = lambda wildcards: "--basename {b}".format(b = wildcards.sample)
    output:
        "bismark_mapped/{sample}_pe.bam", 
	"bismark_mapped/{sample}_PE_report.txt"
    script:
        "wrapper/bismark/align_pe.py"
rule sort_bismark_bam:
    input:
        bismark_mapped_bam = "bismark_mapped/{sample}_pe.bam", 
    params:
    output:
        sorted_bam = "bismark_mapped/{sample}_pe_sorted.bam"
    shell:
        "samtools sort -n {input.bismark_mapped_bam} -o {output.sorted_bam}" 
rule merge_sorted_bam: 
    input:
        bam_file_list = lambda wildcards: config["bam_merge_config"][wildcards.cell_type]
    params:
    output:
        merged_bam = "bismark_mapped/merged_{cell_type}.bam"
    shell:
        "sh scripts/merge_bam.sh \"{input.bam_file_list}\" {output.merged_bam} {wildcards.cell_type}"
    
rule read_to_sam_flag_dict:
    input:
        bam_file = "bismark_mapped/{bam}.bam"
    params:
    output:
        dict_file = "bismark_mapped/alignment_samflag_{bam}.pkl"
        
    shell: 
        "samtools view {input.bam_file} |"
        " awk '{{print $1\"\t\"$2}}' | "
        " python scripts/create_sam_flag_dict.py {output.dict_file}"
