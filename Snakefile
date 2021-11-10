import pandas as pd
import numpy as np
import yaml
from snakemake.utils import validate, min_version
import snakemake


#to run when you add new samples: snakemake -n -R `snakemake --list-input-changes`

##### set minimum snakemake version #####
min_version("5.4.3")

##### Singularity image path
singularity: "/shares/CIBIO-Storage/GROUPS/sharedLC/Davide/containers/tes-analyser-cont.sif"

##### config file
configfile: "configuration/config.yaml"

#read samples info
SAMPLES=pd.read_csv(config['samples'], sep = "\t").set_index("sample", drop=False).sort_index()
#read fastq files
FASTQ_FILES=pd.read_csv(config["fastq_files"], dtype=str, sep ="\t").set_index(["sample", "lane"], drop=False).sort_index()

SAMPLES = set(FASTQ_FILES["sample"])

# print(FASTQ_FILES)
# for row in FASTQ_FILES.itertuples():
#     print(row.sample)
#     print(row.lane)

#define the outputs 
Gene_expression = "results/expression_tabs/Gene_expression_counts.txt"
TEs_expression="results/expression_tabs/TEs_expression_counts.txt"
Gene_TEs_expression="results/expression_tabs/Gene_TEs_expression_counts.txt"
#bigwigs=expand("results/bigwigs/{sample}.bw", sample = SAMPLES)
#deseq2= expand("results/deseq2/{contrast}/log2FC{thresh_fc}_pval{thresh_pval}/{contrast}_enrichments_log2FC{thresh_fc}_pval{thresh_pval}.xlsx", contrast = config['diffexp']['contrasts'],
#                                                                                                                                                thresh_fc = config['diffexp']['log2fc'],
#                                                                                                                                                thresh_pval = config['diffexp']['pval'])

rule all:
    input: 
        Gene_expression,
        TEs_expression,
        Gene_TEs_expression,
        expand("results/bigWigs/{sample}.bw", sample = SAMPLES),
        expand("results/deseq2/{contrast}/log2fc{log2FC}_pval{pvalue}/{contrast}_genes_diffexp_log2fc{log2FC}_pval{pvalue}.tsv", contrast = config["diffexp"]["contrasts"],
            log2FC = config["diffexp"]["log2FC"],
            pvalue = config["diffexp"]["pval"]),
        expand("results/deseq2/{contrast}/log2fc{log2FC}_pval{pvalue}/{contrast}_volcano_genes_log2fc{log2FC}_pval{pvalue}.pdf", contrast = config["diffexp"]["contrasts"],
            log2FC = config["diffexp"]["log2FC"],
            pvalue = config["diffexp"]["pval"]),
        expand("results/deseq2/{contrast}/log2fc{log2FC}_pval{pvalue}/{contrast}_enrichments_log2fc{log2FC}_pval{pvalue}.xlsx", contrast = config["diffexp"]["contrasts"],
            log2FC = config["diffexp"]["log2FC"],
            pvalue = config["diffexp"]["pval"]),
        expand("results/deseq2/{contrast}/log2fc{log2FC}_pval{pvalue}/{contrast}_MA_plot_tes_log2fc{log2FC}.pdf", contrast = config["diffexp"]["contrasts"],
            log2FC = config["diffexp"]["log2FC"],
            pvalue = config["diffexp"]["pval"])

#load the rules 
include: "workflow/rules/dag.smk"
include: "workflow/rules/trimming.smk"
include: "workflow/rules/common.smk"
include: "workflow/rules/align.smk"
include: "workflow/rules/filtering_bam.smk"
include: "workflow/rules/diffExpr.smk"

##### handle possible errors, clean temp folders #####
onsuccess:
    # shell("""
    # #rm -r /results/fastq/
    # """)
    print("All done!\n")

onerror:
    print("An error ocurred. Workflow aborted")