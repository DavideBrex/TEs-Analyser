
def return_fastq_files(wildcars):
    return units.loc[(wildcards.sample, wildcards.lane), ["fq1", "fq2"]].dropna()



rule link_to_fastq_se:
    input: 
        return_fastq_files
    output:
        fastq_out = temp("/results/fastq/{sample}-{lane}.fastq.gz") 
    message:
        "Creating link to fastq file {input}"
    shell:
        """
        ln -s {input} {output}
        """

rule link_to_fastq_pe:
    input: r
        return_fastq_files
    output:
        fastq_out_fw = temp("/results/fastq/{sample}-{lane}.1.fastq.gz") 
        fastq_out_rv = temp("/results/fastq/{sample}-{lane}.2.fastq.gz") 
    message:
        "Creating link to fastq file {input}"
    shell:
        """
        ln -s {input[0]} {output.fastq_out_fw}
        ln -s {input[1]} {output.fastq_out_rv}
        """


rule merge_fastq_pe:
    input:
        fw = lambda w: expand("results/fastq/{lane.sample}-{lane.lane}.1.fastq.gz", lane=units.loc[w.sample].itertuples()),
        rv = lambda w: expand("results/fastq/{lane.sample}-{lane.lane}.2.fastq.gz", lane=units.loc[w.sample].itertuples())
    output:
        fastq1 = temp("results/fastq/{{sample}}.1.fastq.gz"),
        fastq2 = temp("results/fastq/{{sample}}.2.fastq.gz")
    log:
        "results/logs/fastp/{sample}.log"
    message:
        "Merging fastq files from {input}"
    shell:
        """
        cat {input.fw} > {output.fastq1}
        cat {input.rv} > {output.fastq2}
        """


rule merge_fastq_se:
    input:
        lambda w: expand("results/fastq/{lane.sample}-{lane.lane}.fastq.gz", lane=units.loc[w.sample].itertuples()),
    output:
        temp("results/fastq/{{sample}}.se.fastq.gz")
    log:
        "results/logs/fastp/{sample}.log"
    message:
        "Merging fastq files from {input}"
    shell:
        """
        cat {input} > {output}
        """



rule fastp_pe:
    input:
        fw = "results/fastq/{{sample}}.1.fastq.gz",
        rv = "results/fastq/{{sample}}.2.fastq.gz"
    output:
        fastq1 = temp("results/fastq/trimmed/{{sample}}.1.fastq.gz"),
        fastq2 = temp("results/fastq/trimmed/{{sample}}.2.fastq.gz")
    log:
        "logs/fastp/{sample}.log"
    threads:
        config["tools_cpu"]["fastp_pe"]
    params:
        fastp_params = config["params"]["fastp"]["pe"],
    message:
        "Processing fastq files from {input}"
    shadow:
        "minimal"
    shell:
        """
        fastp -i {input.fw} \
        -I {input.rv} \
        -o {output.fastq1} \
        -O {output.fastq2} \
        -w {threads} \
        {params.fastp_params} 2> {log}
        """


rule fastp_se:
    input:
        "results/fastq/{{sample}}.se.fastq.gz"
    output:
        temp("results/fastq/trimmed/{{sample}}.se.fastq.gz")
    log:
        "logs/fastp/{sample}.log"
    threads:
        config["tools_cpu"]["fastp_pe"]
    params:
        fastp_params = config["params"]["fastp"]["se"],
    message:
        "Processing fastq files from {input}"
    shadow:
        "minimal"
    shell:
        """
        fastp -i {input} \
        -o {output} \
        -w {threads} {params.fastp_params} 2> {log}
        """