
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
    log:
        "results/logs/bam_tagging/{sample}.log"
    message:
        "Tagging the bam file {input}..."
    shell:
        """
        java -Xmx4g -jar {params.dropseq} \
        TagReadWithGeneFunction I={input} \
        O={output} \
        ANNOTATIONS_FILE={params.gtf_file} \
        USE_STRAND_INFO={params.dropseq_use_strand} 2> {log}
        """

# select multimappers reads and intronic + intergenics reads 
rule filter_bam_file:
    input: 
        rules.tag_bam_file.output
    output: 
        temp("results/allbams/{sample}.filtered.bam")
    params:
        filtering_json = config["params"]["bam_filtering_json"]
    log:
        "results/logs/bam_filtering/{sample}.log"
    message:
        "Selecting multimapping, intronic or intergenic reads..."
    shell:
        """
        bamtools filter -in {input} \
        -script {params.filtering_json} \
        -out {output} 2> {log}
        """


#convert bam files to fastq files
# we add  awk to to the following:
# since in the fastq file we would have multiple times the same read id (and mate), we keep 
# only one read id (2 actually, because two mates) and remove all the other duplicated
# in this way the resulting fastq is lighter and the re-alignemnts is faster. 

rule bam_to_fastq:
    input:
        rules.filter_bam_file.output
    output:
        "results/filtered_fastq/{sample}.filtered.fastq.gz"
    threads:
        config["tools_cpu"]["samtools_sort"]
    params:
        "results/filtered_fastq/{sample}"
    log:
        "results/logs/bam_to_fastq/{sample}.log"
    message:
        "Converting the filtered bam file {input} to fastq..."
    shell:
        """
        samtools sort -n -T {params}.tmp -@ {threads} {input} -o {params}.bam;

        samtools view -H {params}.bam > {params}.tmp.bam \
        && samtools view {params}.bam | awk '!_[$1,$10]++' >> {params}.tmp.bam \
        && samtools view -b {params}.tmp.bam > {params}.bam && rm {params}.tmp.bam;     
    

        bedtools bamtofastq -i {params}.bam -fq /dev/stdout \
        | gzip -c > {output} && rm {params}.bam 2> {log}
        """ 
