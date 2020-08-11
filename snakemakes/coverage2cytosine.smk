rule coverage2cytosine:
    input:
        "bismark_methylation_extracted/{sample}_pe.bismark.cov.gz"
    params:
        important_params = "--gc --nome-seq --zero_based --gzip",
        genome_folder = "/beevol/home/satyanarr/workplace/data/ucsc/dm/dm3"
    output:
        "coverage_to_cytosine/{sample}.NOMe.CpG.cov.gz", 
        "coverage_to_cytosine/{sample}.NOMe.GpC.cov.gz", 
        "coverage_to_cytosine/{sample}.NOMe.CpG_report.txt.gz", 
        "coverage_to_cytosine/{sample}.NOMe.GpC_report.txt.gz" 
    shell:
        "coverage2cytosine"
        " --dir coverage_to_cytosine --genome_folder {params.genome_folder}"
        " {params.important_params} -o {wildcards.sample} {input}"
