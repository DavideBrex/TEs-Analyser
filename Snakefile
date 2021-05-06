import pandas as pd
import numpy as np
import yaml
from snakemake.utils import validate, min_version
import snakemake

##### set minimum snakemake version #####
min_version("5.4.3")

##### Singularity image path
#singularity: "/hpcnfs/data/DP/Singularity/dfernandezperez-bioinformatics-singularity-master-chipseq.simg"

##### config file
configfile: "configuration/config.yaml"

#read samples info
SAMPLES=pd.read_csv(config['samples'], sep = "\t").set_index("sample", drop=False).sort_index()
#read fastq files
FASTQ_FILES=pd.read_csv(config["fastq_files"], dtype=str, sep ="\t").set_index(["sample", "lane"], drop=False).sort_index()

SAMPLES = set(FASTQ_FILES["sample"])

print(FASTQ_FILES)

#define the outputs 
Gene_expression = "results/expression_tabs/Gene_expression_counts.txt"
TEs_expression="results/expression_tabs/TEs_expression.txt"
bigwigs=expand("results/bigwigs/{sample}.bw", sample = SAMPLES)
deseq2= expand("results/deseq2/{contrast}/log2FC{thresh_fc}_pval{thresh_pval}/{contrast}_enrichments_log2FC{thresh_fc}_pval{thresh_pval}.xlsx", contrast = config['diffexp']['contrasts'],
                                                                                                                                                thresh_fc = config['diffexp']['log2fc'],
                                                                                                                                                thresh_pval = config['diffexp']['pval'])

rule all:
    input: Gene_expression + TEs_expression  + deseq2


#load the rules 
include: "rules/common.smk"
include: "rules/align.smk"
include: "rules/trimming.smk"



##### handle possible errors, clean temp folders #####
onsuccess:
    shell("""
    rm -r fastq/
    """)
    print("All done!\n")

onerror:
    print("An error ocurred. Workflow aborted")