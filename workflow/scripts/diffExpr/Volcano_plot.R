#set log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(ggplot2))

#read tables
table_genes <- read.delim(snakemake@input[["table_genes_ann"]])
table_tes <- read.delim(snakemake@input[["table_tes_ann"]])
table_genes_tes <- read.delim(snakemake@input[["table_genes_tes_ann"]])

#read parameters
log2FC   <- as.numeric(snakemake@params[["log2fc"]])
pval      <- as.numeric(snakemake@params[["pval"]])
plot_title <- gsub("-", " ", snakemake@params[["contrast"]])

#plotting function
plot_volcano_fun <- function(df, output_name, title, pval, log2FC){
  
    # Credits for the code: Daniel Fernandez Perez
    volcano_plot <- ggplot(data = df, aes(x=log2FoldChange, y=-log10(padj), colour=DEG) ) +
                        geom_point(alpha=0.7, size=2) +

                        annotate("text", label = sum(df$DEG == "Upregulated"), color = "red", y = 0, x = 9, 
                                    vjust="inward",hjust="inward", size = 7) +
                        annotate("text", label = sum(df$DEG == "Downregulated"), color = "darkgreen", y = 0, x = -8,
                                    vjust="inward",hjust="inward", size = 7) +

                        theme_classic() +
                        theme(legend.title = element_blank()) +
                        theme(legend.position = "top") +

                        ggtitle(title) +
                        theme(plot.title = element_text(lineheight=.8, face="bold", hjust = .5)) +

                        xlim(c(-8,8)) + ylim(c(0,30)) + 

                        geom_hline(yintercept = -log10(pval), linetype = 2) +
                        geom_vline(xintercept = c(-log2FC, log2FC), linetype = 2) +

                        xlab("log2 fold change") + ylab("-log10 p-value") +
                        scale_colour_manual(values=c("darkgreen", "gray", "red") )
    

    pdf(snakemake@output[[output_name]])
    volcano_plot
    dev.off()
}

#generate the volcano plots

volcano_genes <- plot_volcano_fun(table_genes,
                                  output_name = "volcano_genes",
                                  title = paste0("Genes - ", plot_title), pval, log2FC)
volcano_tes <- plot_volcano_fun(table_tes, 
                                output_name = "volcano_tes",
                                title = paste0("TEs - ", plot_title), pval, log2FC)
volcano_genes_tes <- plot_volcano_fun(table_genes_tes, 
                                      output_name = "volcano_genes_tes",
                                      title = paste0("Genes and TEs - ", plot_title), pval, log2FC)


