



rule align_first_pass:
    input: get_fq
    output: 
        bam = "results/alignments/{sample}.bam",
        index = "results/alignments/{sample}.bam.bai"
    threads:
        config["threads"]["align"]
    params:
        genome_index = config["ref"]["genome_index"]
        star_par=config["params"]["star_first_pass"]
    message:
        "First-pass alignment on the reference genome for: {input}"
    log:
        align = "results/logs/alignments/{sample}.log"
    shell:
        """
        STAR --genomeDir {params.genome_index} --readFilesIn {input} \
        --runThreadN {threads} {params.star_par}
        """