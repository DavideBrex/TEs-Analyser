#set log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

message("Loading libraries...\n")
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(tible))

#read tables
table_genes <- read.delim(snakemake@input[["table_genes"]])
table_tes <- read.delim(snakemake@input[["table_tes"]])
table_genes_tes <- read.delim(snakemake@input[["table_genes_tes"]])

#read parameters
log2FC   <- as.numeric(snakemake@params[["log2fc"]])
pval      <- as.numeric(snakemake@params[["pval"]])

# add label to Up and Down regulated genes
Annot_DE <- function(df, log2FC = 2, padjust = 0.05) {
  df$DEG <- "NotDE"
  df$DEG[which(df$log2FoldChange >= log2FC & df$padj <= padjust)] <- "Upregulated" #Annotate UPregulated genes
  df$DEG[which(df$log2FoldChange <= -log2FC & df$padj <= padjust)] <- "Downregulated" #Annotate DOWNregulated genes
  df <- df %>% dplyr::select(-lfcSE, -pvalue)
  return(df)
}

table_genes_ann <- Annot_DE(table_genes, log2FC, pval)
table_tes_ann <- Annot_DE(table_tes, log2FC, pval)
table_genes_tes_ann <- Annot_DE(table_genes_tes, log2FC, pval)


#store output
write.table(table_genes_ann, file = snakemake@output[["table_genes_ann"]], sep = "\t",quote = F, row.names = F)
write.table(table_tes_ann, file = snakemake@output[["table_tes_ann"]], sep = "\t",quote = F, row.names = F)
write.table(table_genes_tes_ann, file = snakemake@output[["table_genes_tes_ann"]], sep = "\t",quote = F, row.names = F)




