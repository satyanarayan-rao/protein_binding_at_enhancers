rule extract_methylation_info_from_bam:
    input:
        "bismark_mapped/{sample}_pe.bam"
    params:
        genome_folder = "/beevol/home/satyanarr/workplace/data/ucsc/dm/dm3",
        extraction = "--parallel 4 --cytosine_report",
        extra = "--buffer_size 100% --CX --bedGraph"
         
    output:
        "bismark_methylation_extracted/CpG_context_{sample}_pe.txt.gz",
        "bismark_methylation_extracted/CHG_context_{sample}_pe.txt.gz",
        "bismark_methylation_extracted/CHH_context_{sample}_pe.txt.gz",
        "bismark_methylation_extracted/{sample}_pe.bismark.cov.gz",
        "bismark_methylation_extracted/{sample}_pe.bedGraph.gz",
        "bismark_methylation_extracted/{sample}_pe.CX_report.txt.gz"
    shell:
        "bismark_methylation_extractor -p --comprehensive"
        " --genome_folder {params.genome_folder} --gzip -o"
        " bismark_methylation_extracted {params.extraction} {params.extra} {input}"
