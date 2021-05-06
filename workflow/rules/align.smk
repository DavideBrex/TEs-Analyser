



rule align_first_pass:
    input: get_fq
    output: 
        bam = "results/alignments/{sample}.bam",
        index = "results/alignments/{sample}.bam.bai"
    threads:
        config["tools_cpu"]["STAR_first_pass"]
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
        --runThreadN {threads} --quantMode GeneCounts --readFilesCommand zcat \   
        --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS NM MD  {params.star_par}
        """


rule align_second_pass:
    input: 