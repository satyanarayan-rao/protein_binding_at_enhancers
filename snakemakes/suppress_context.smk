rule suppress_context:
    input:
        bam_file = "bismark_mapped/merged_S2.bam"
    output:
        out_bam = "bismark_mapped/suppressed_merged_S2.bam"
    shell:
        "sh scripts/suppress_context.sh {input.bam_file} {output.out_bam}"
    
