

#perform differential expression analysis
rule create_deseq2_ojects:
    input: 
        tes_expr = rules.TEs_counting.output.tes_expression_counts,
        genes_expr = rules.merge_count_tables.output,
        genes_tes_expr = rules.TEs_counting.output.gene_tes_expression_counts
    output: 
        rds_genes = "results/deseq2/deseq2_object_genes.rds",
        rds_tes = "results/deseq2/deseq2_object_TEs.rds",
        rds_genes_tes = "results/deseq2/deseq2_object_genes_TEs.rds"
    params:
        sample_file = config["samples"],
        pval_thres = config["diffexp"]["pval"],
        log2fc = config["diffexp"]["log2fc"]
    log:
        "results/logs/deseq2/deseq2_objects.log"
    script:
        "../scripts/diffExpr/Deseq2_objects.R" 



rule differential_expresssion:
    input:
        rds_genes = rules.create_deseq2_ojects.output.rds_genes,
        rds_tes = rules.create_deseq2_ojects.output.rds_tes,
        rds_genes_tes = rules.create_deseq2_ojects.output.rds_genes_tes
    output:
        table_genes = "results/deseq2/{contrast}/{contrast}_genes_diffexp.tsv",
        table_tes = "results/deseq2/{contrast}/{contrast}_tes_diffexp.tsv",
        table_genes_tes = "results/deseq2/{contrast}/{contrast}_genes_tes_diffexp.tsv"
    params:
        contrast = lambda w: config["diffexp"]["contrasts"][w.contrast]
    log:
        "results/logs/deseq2/{contrast}.diffexp.log"
    script:
        "../scripts/diffExpr/Deseq2_DiffExpr.R"


rule filtering_de_elements:
    input: 
        table_genes = rule.differential_expresssion.output.table_genes,
        table_tes = rule.differential_expresssion.output.table_tes,
        table_genes_tes = rule.differential_expresssion.output.table_genes_tes
    output: 
        table_genes_ann = "results/deseq2/{contrast}/log2fc{log2FC}_pval{pvalue}/{contrast}_genes_diffexp_log2fc{log2FC}_pval{pvalue}.tsv",
        table_tes_ann = "results/deseq2/{contrast}/log2fc{log2FC}_pval{pvalue}/{contrast}_tes_diffexp_log2fc{log2FC}_pval{pvalue}.tsv",
        table_genes_tes_ann = "results/deseq2/{contrast}/log2fc{log2FC}_pval{pvalue}/{contrast}_genes_tes_diffexp_log2fc{log2FC}_pval{pvalue}.tsv"
    params:
        pval      = lambda w: w.pvalue,
        log2fc    = lambda w: w.log2FC
    log:
        "results/logs/deseq2/{contrast}.{pvalue}.{log2fc}.filter_de.log"
    script:
        "../scripts/diffExpr/Filter_DE_elements.R"



rule volcano:
    input:
        table_genes_ann = rules.filtering_de_elements.output.table_genes_ann,
        table_tes_ann = rules.filtering_de_elements.output.table_tes_ann,
        table_genes_tes_ann = rules.filtering_de_elements.output.table_genes_tes_ann
    output:
        volcano_genes = "results/deseq2/{contrast}/log2fc{log2FC}_pval{pvalue}/{contrast}_volcano_genes_log2fc{log2FC}_pval{pvalue}.pdf",
        volcano_tes = "results/deseq2/{contrast}/log2fc{log2FC}_pval{pvalue}/{contrast}_volcano_tes_log2fc{log2FC}_pval{pvalue}.pdf",
        volcano_genes_tes = "results/deseq2/{contrast}/log2fc{log2FC}_pval{pvalue}/{contrast}_volcano_genes_tes_log2fc{log2FC}_pval{pvalue}.pdf"
    params:
        pval     = lambda w: w.pvalue,
        log2fc   = lambda w: w.log2fc,
        contrast = lambda w: w.contrast
    log:
        "results/logs/deseq2/{contrast}.pval_{pvalue}.log2fc_{log2fc}.volcano.log"
    script:
        "../scripts/diffExpr/Volcano_plot.R"  

rule enrichments:
    input:
        rules.filtering_de_elements.output.table_genes_ann
    output:
        "results/deseq2/{contrast}/{contrast}_enrichments.xls"
    params:
        genome = config["ref"]["genome"]
    log:
        "results/logs/deseq2/{contrast}.GSEA.log"
    script:
        "../scripts/diffExpr/GSEA_analysis.R"   

