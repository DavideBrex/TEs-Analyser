
# use dropseq to annotate the reads in the bam files
rule tag_bam_file:
    input: 
        "results/allbams/{sample}.bam"
    output: 
        temp("results/allbams/{sample}.tagged.bam")
    params:
        gtf_file = config["ref"]["gtf_file"],
        dropseq = config["dropseq_jar"],
        dropseq_use_strand = config["params"]["drop_seq_strand"]
    message:
        "Tagging the bam file {input}..."
    shell:
        """
        java -Xmx4g -jar {params.dropseq} \
        TagReadWithGeneFunction I={input} \
        O={output} \
        ANNOTATIONS_FILE={params.gtf_file} \
        USE_STRAND_INFO={params.dropseq_use_strand}
        """

# select multimappers reads and intronic + intergenics reads 
rule filter_bam_file:
    input: 
        rules.tag_bam_file.output
    output: 
        temp("results/allbams/{sample}.filtered.bam")
    params:
        filtering_json = config["params"]["bam_filtering_json"]
    message:
        "Selecting multimapping, intronic or intergenic reads..."
    shell:
        """
        bamtools filter -in {input} \
        -script {params.filtering_json} \
        -out {output}
        """


#convert bam files to fastq files
rule bam_to_fastq:
    input:
        rules.filter_bam_file.output
    output:
        "results/filtered_fastq/{sample}.filtered.fastq.gz"
    threads:
        config["tools_cpu"]["samtools_sort"]
    params:
        "results/filtered_fastq/{sample}"
    message:
        "Converting the filtered bam file {input} to fastq..."
    shell:
        """
        samtools sort -n -T {params}.tmp -@ {threads} {input} \
        | bedtools bamtofastq -i stdin -fq /dev/stdout \
        | gzip -c > {output} 
        """