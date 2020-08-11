rule trim_galore_pe:
    input:
        #["raw_data/{sample}_1.fastq.gz", "raw_data/{sample}_2.fastq.gz"]
        lambda wildcards: get_raw_fasta(wildcards) 
    params:
        #extra = "-q 30 --illumina --gzip --trim-n -o trim_galore_out --cores 3"
        extra = lambda wildcards: "-q 30 --illumina --gzip --trim-n -o trim_galore_out --cores 3 --no_report_file --basename {s}".format (s = wildcards.sample)
    output:
        #"trimmed/{sample}_1_val_1.fq.gz", 
        #"trimmed/{sample}_1.fastq.gz_trimming_report.txt",
        #"trimmed/{sample}_2_val_2.fq.gz",
        #"trimmed/{sample}_2.fastq.gz_trimming_report.txt"
        "trimmed/{sample}_R1_val_1.fq.gz", 
        "trimmed/{sample}_R2_val_2.fq.gz",

    log:
        "logs/trim_galore/{sample}.log"
    script:
        #"0.42.0/bio/trim_galore/pe"
        "wrapper/trim_galore/pe.py"
