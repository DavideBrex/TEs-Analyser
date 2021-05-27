message("Loading libraries...\n")
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(IRanges))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(DESeq2))

#set log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


#read the sample information
metadata <- read.table(snakemake@params[["sample_file"]], header = T)

#read the count tables
tes_counts <- read.table(snakemake@input[["tes_expr"]], header = T, sep = "\t", check.names = F)
genes_counts <- read.table(snakemake@input[["genes_expr"]], header = T, sep = "\t", check.names = F)
tes_genes_counts <- read.table(snakemake@input[["genes_tes_expr"]], header = T, sep = "\t", check.names = F)

tes_counts <- column_to_rownames(tes_counts, var= "TE_id")
genes_counts <- column_to_rownames(genes_counts, var= "Gene_id")
tes_genes_counts <- column_to_rownames(tes_genes_counts, var= "Gene_TE_id")

# sort the columns in the count tables so that are in the same order as in metadata
tes_counts <- tes_counts[, as.character(metadata$sample)]
stopifnot(identical(colnames(tes_counts), as.character(metadata$sample)))

genes_counts <- genes_counts[, as.character(metadata$sample)]
stopifnot(identical(colnames(genes_counts), as.character(metadata$sample)))

tes_genes_counts <- tes_genes_counts[, as.character(metadata$sample)]
stopifnot(identical(colnames(tes_genes_counts), as.character(metadata$sample)))

#create deseq2 objects

tes_counts_dds <- DESeqDataSetFromMatrix(tes_counts, metadata, design = ~ condition)
tes_dds <- DESeq(tes_counts_dds)

genes_counts_dds <- DESeqDataSetFromMatrix(genes_counts, metadata, design = ~ condition)
genes_dds <- DESeq(genes_counts_dds)

genes_tes_counts_dds <- DESeqDataSetFromMatrix(tes_genes_counts, metadata, design = ~ condition)
genes_tes_dds <- DESeq(genes_tes_counts_dds)

# store the objects

saveRDS(tes_dds, file = snakemake@output[["rds_tes"]])
saveRDS(genes_dds, file = snakemake@output[["rds_genes"]])
saveRDS(genes_tes_dds, file = snakemake@output[["rds_genes_tes"]])
