


def is_single_end(sample):
    #check whether the sample is single end
    return pd.isnull(FASTQ_FILES.loc[(sample), "fq2"][0])

def get_fq(wildcards):
    """
    Return the raw fastq files
    """
    if config["trimming"]:
        if not is_single_end(**wildcards):
            #paired-end sample
            return expand("results/fastq/trimmed/{sample}.{group}.fastq.gz", group=[1, 2], **wildcards)
        #single end sample
        return "results/fastq/trimmed/{sample}.se.fastq.gz", **wildcards)
    else:
        # no trimming, use raw reads
        if not is_single_end(**wildcards):
            #paired-end sample
            return expand("results/fastq/{sample}.{group}.fastq.gz", group=[1, 2], **wildcards)
        #single end sample  
        return "results/fastq/{sample}.se.fastq.gz", **wildcards)


def get_modified_fq(wildcards):
    """
    Return the newly created fastq files (were created from the filtered bam files)
    """
    return "results/filtered_fastq/{sample}.fastq.gz", **wildcards)

