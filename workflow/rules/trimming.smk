
def return_fastq_files(wildcards):
    return FASTQ_FILES.loc[(wildcards.sample, wildcards.lane), ["fq1", "fq2"]].dropna()

def return_lanes_forward(wildcards):
    return expand("results/fastq/{lane.sample}-{lane.lane}.1.fastq.gz", lane=FASTQ_FILES.loc[wildcards.sample].itertuples())

def return_lanes_reverse(wildcards):
    return expand("results/fastq/{lane.sample}-{lane.lane}.2.fastq.gz", lane=FASTQ_FILES.loc[wildcards.sample].itertuples())

#create a symbolic link to fastq files (single end)
rule copy_fastq_se:
    input: 
        return_fastq_files
    output:
        fastq_out = temp("results/fastq/{sample}-{lane}.fastq.gz") 
    message:
        "Creating link to fastq file {input}"
    shell:
        """
        ln -s {input} {output}
        """

#create a symbolic link to fastq files (paired end)
rule copy_fastq_pe:
    input: 
        return_fastq_files
    output:
        fastq_out_fw = temp("results/fastq/{sample}-{lane}.1.fastq.gz"),
        fastq_out_rv = temp("results/fastq/{sample}-{lane}.2.fastq.gz") 
    message:
        "Creating link to fastq file {input}"
    shell:
        """
        ln -s {input[0]} {output.fastq_out_fw}
        ln -s {input[1]} {output.fastq_out_rv}
        """

#merge the multiple lines if present (for paired end)
rule merge_fastq_pe:
    input:
        fw = return_lanes_forward,
        rv = return_lanes_reverse
    output:
        fastq1 = temp("results/fastq/{sample}.1.fastq.gz"),
        fastq2 = temp("results/fastq/{sample}.2.fastq.gz")
    log:
        "results/logs/fastp/{sample}.log"
    message:
        "Merging fastq files from {input}"
    shell:
        """
        cat {input.fw} > {output.fastq1}
        cat {input.rv} > {output.fastq2}
        """

#merge the multiple lines if present (for single end)
rule merge_fastq_se:
    input:
        lambda w: expand("results/fastq/{lane.sample}-{lane.lane}.fastq.gz", lane=FASTQ_FILES.loc[w.sample].itertuples()),
    output:
        temp("results/fastq/{sample}.se.fastq.gz")
    log:
        "results/logs/fastp/{sample}.log"
    message:
        "Merging fastq files from {input}"
    shell:
        """
        cat {input} > {output}
        """


# run fastp for paired end samples
rule fastp_pe:
    input:
        fw = "results/fastq/{sample}.1.fastq.gz",
        rv = "results/fastq/{sample}.2.fastq.gz"
    output:
        fastq1 = temp("results/fastq/trimmed/{sample}.1.fastq.gz"),
        fastq2 = temp("results/fastq/trimmed/{sample}.2.fastq.gz")
    log:
        "logs/fastp/{sample}.log"
    threads:
        config["tools_cpu"]["fastp_pe"]
    params:
        fastp_params = config["params"]["fastp"]["pe"],
    message:
        "Processing fastq files from {input}"
    shell:
        """
        fastp -i {input.fw} \
        -I {input.rv} \
        -o {output.fastq1} \
        -O {output.fastq2} \
        -w {threads} \
        {params.fastp_params} 2> {log}
        """

# run fastp for single end samples
rule fastp_se:
    input:
        "results/fastq/{sample}.se.fastq.gz"
    output:
        temp("results/fastq/trimmed/{sample}.se.fastq.gz")
    log:
        "logs/fastp/{sample}.log"
    threads:
        config["tools_cpu"]["fastp_pe"]
    params:
        fastp_params = config["params"]["fastp"]["se"],
    message:
        "Processing fastq files from {input}"
    shell:
        """
        fastp -i {input} \
        -o {output} \
        -w {threads} {params.fastp_params} 2> {log}
        """