

#perform differential expression analysis
rule create_deseq2_ojects:
    input: 
        tes_expr = rules.TEs_counting.output.tes_expression_counts
        genes_expr = rules.merge_count_tables.output
        genes_tes_expr = rules.TEs_counting.output.gene_tes_expression_counts
    output: 
        rds_genes = "results/deseq2/deseq2_object_genes.rds"
        rds_tes = "results/deseq2/deseq2_object_TEs.rds"
        rds_genes_tes = "results/deseq2/deseq2_object_genes_TEs.rds"
    params:
        sample_file = config["samples"]
        pval_thres = config["diffexp"]["pval"]
        log2fc = config["diffexp"]["log2fc"]
    log:
        "results/logs/deseq2/deseq2_objects.log"
    script:
        "../scripts/diffExpr/Deseq2_objects.R" 



rule differential_expresssion:
    input:
        rds_genes = rules.create_deseq2_ojects.output.rds_genes
        rds_tes = rules.create_deseq2_ojects.output.rds_tes
        rds_genes_tes = rules.create_deseq2_ojects.output.rds_genes_tes
    output:
        table_genes = "results/deseq2/{contrast}/{contrast}_genes_diffexp.tsv"
        table_tes = "results/deseq2/{contrast}/{contrast}_tes_diffexp.tsv"
        table_genes_tes = "results/deseq2/{contrast}/{contrast}_genes_tes_diffexp.tsv"
    params:
        contrast = lambda w: config["diffexp"]["contrasts"][w.contrast]
    log:
        "results/logs/deseq2/deseq2_diffexp.log"
    script:
        "../scripts/diffExpr/Deseq2_DiffExpr.R"


rule filtering_de_elements:
    input: rule.create_deseq2_ojects.output
    output: 
    run: 