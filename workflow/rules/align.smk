
rule align_first_pass:
    input: get_fq
    output: 
        bam = "results/alignments/{sample}/{sample}.bam",
        index = "results/alignments/{sample}/{sample}.bam.bai"
        log   = "results/alignments/{sample}/Log.final.out"
        counts = "results/alignments/{sample}/ReadsPerGene.out.tab"
    threads: config["tools_cpu"]["STAR_first_pass"]
    params:
        genome_index = config["ref"]["genome_index"]
        star_par=config["params"]["star_first_pass"]
    message:
        "First-pass alignment on the reference genome for: {input}"
    log:
        align = "results/logs/alignments/{sample}.log"
    shell:
        """
        STAR --genomeDir {params.genome_index} \
        --readFilesIn {input} \
        --runThreadN {threads} \
        --quantMode GeneCounts \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS NM MD  \
        {params.star_par} 2> {log.align} 
        samtools index {output.bam} 2>> {log.align}
        """


rule align_second_pass:
    input: get_modified_fq
    output: 
        bam = "results/alignments/{sample}-filtered/{sample}.bam",
        index = "results/alignments/{sample}-filtered/{sample}.bam.bai"
        log   = "results/02alignments/{sample}-filtered/{sample}/Log.final.out"
    threads: config["tools_cpu"]["STAR_second_pass"]
    params:
        pseudo_genome_index = config["ref"]["pseudo_genome_index"]
        star_par=config["params"]["star_second_pass"]
    message:
        "First-pass alignment on the reference genome for: {input}"
    log:
        align = "results/logs/alignments/{sample}.log"
    shell:
        """
        STAR --genomeDir {params.genome_index} \
        --readFilesIn {input} \
        --runThreadN {threads} \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS NM MD  \
        --alignIntronMax 1 \
        {params.star_par} 2> {log.align} 
        samtools index {output.bam} 2>> {log.align}
        """

#copy the bam files to the folder allbams. remove it when all the pipeline is done
rule move_bams:
    input: rules.align_first_pass.output.bam
    output: bam = temp("results/allbams/{sample}.bam")
    shell:
        """
        cp {input} {output}
        """


rule merge_count_tables:
    input: 
        
    output: 
    run: 