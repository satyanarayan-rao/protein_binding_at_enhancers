rule bismark2report: 
    input: 
        "bismark_mapped/{sample}_PE_report.txt"
    params: 
        title = lambda wildcards: wildcards.sample,
        basename = lambda wildcards: "{d}/{b}".format(d="bismark_summary", b = wildcards.sample)
    output: 
        "bismark_align_report/{sample}.html"
    shell:
        "bismark2report --output {output} --alignment_report {input}" 
