

rule align_first_pass:
    input: get_fq
    output: 
        bam = "results/alignments/{sample}/{sample}.bam",
        index = "results/alignments/{sample}/{sample}.bam.bai",
        log   = "results/alignments/{sample}/Log.final.out",
        counts = "results/alignments/{sample}/ReadsPerGene.out.tab"
    threads: config["tools_cpu"]["STAR_first_pass"]
    params:
        genome_index = config["ref"]["genome_index"],
        star_par = config["params"]["star_first_pass"],
	out_dir = "results/alignments/{sample}/",
	samtools_cpus = config["tools_cpu"]["samtools_sort"]
    message:
        "First-pass alignment on the reference genome for: {input}"
    log:
        align = "results/logs/alignments/{sample}.firstpass.log"
    shell:
        """
        STAR --genomeDir {params.genome_index} \
        --readFilesIn {input} \
        --runThreadN {threads} \
        --quantMode GeneCounts \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS NM MD  \
	--outFileNamePrefix {params.out_dir} \
        {params.star_par} 2> {log.align}
	mv {params.out_dir}Aligned.sortedByCoord.out.bam {output.bam}
        samtools index {output.bam} 2>> {log.align}
        """


rule align_second_pass:
    input: get_modified_fq
    output: 
        bam = "results/alignments/{sample}-filtered/{sample}.bam",
        index = "results/alignments/{sample}-filtered/{sample}.bam.bai",
        log   = "results/alignments/{sample}-filtered/Log.final.out"
    threads: config["tools_cpu"]["STAR_second_pass"]
    params:
        pseudo_genome_index = config["ref"]["pseudo_genome_index"],
        star_par = config["params"]["star_second_pass"],
	out_dir = "results/alignments/{sample}-filtered/"
    message:
        "Second-pass alignment on the reference genome for: {input}"
    log:
        align = "results/logs/alignments/{sample}.secondpass.log"
    shell:
        """
        STAR --genomeDir {params.pseudo_genome_index} \
        --readFilesIn {input} \
        --runThreadN {threads} \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS NM MD  \
        --alignIntronMax 1 \
	--outFileNamePrefix {params.out_dir} \
        {params.star_par} 2> {log.align}
	mv {params.out_dir}Aligned.sortedByCoord.out.bam {output.bam} 
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


#combine the single sample gene count tables
rule merge_count_tables:
    input: 
        expand("results/alignments/{sample}/ReadsPerGene.out.tab", sample = SAMPLES)
    output: 
        "results/expression_tabs/Gene_expression_counts.txt"
    params:
        col_to_pick=config["params"]["htseq_count_column"]
    script:
        "../scripts/merge_count_tables.R"


rule TEs_counting:
    input: 
        bam_second_pass = expand("results/alignments/{sample}-filtered/{sample}.bam", sample = SAMPLES),
        gene_count_expression = "results/expression_tabs/Gene_expression_counts.txt"
    output:
        tes_expression_counts = "results/expression_tabs/TEs_expression_counts.txt",
        gene_tes_expression_counts = "results/expression_tabs/Gene_TEs_expression_counts.txt"
    params:
        gtf_TEs = config["ref"]["gtf_file_TEs"]
    script:
        "../scripts/quantify_TEs_expression.R" 
