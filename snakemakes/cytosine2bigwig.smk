rule cytosine2bigwig:
    input: 
        "coverage_to_cytosine/{sample}.cov.gz"
    params:
        chrom_size = "/beevol/home/satyanarr/data/ucsc/dm/dm3/dm3.chrom.sizes"
    output:
        "bigwigs/{sample}.bigwig"
    shell:
        "sh scripts/cytosine2bigwig.sh {input} {output} {params.chrom_size}"
