#set log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(ggrepel))
suppressPackageStartupMessages(require(tidyverse))


#read tables
table_genes <- read.delim(snakemake@input[["table_genes_ann"]])
table_tes <- read.delim(snakemake@input[["table_tes_ann"]])
table_genes_tes <- read.delim(snakemake@input[["table_genes_tes_ann"]])

#read parameters
log2FC   <- as.numeric(snakemake@params[["log2FC"]])
pval      <- as.numeric(snakemake@params[["pvalue"]])
plot_title <- gsub("-", " ", snakemake@params[["contrast"]])

#plotting function
plot_volcano_fun <- function(df, output_name, title, pval, log2FC){
  
    # Credits for the code: Daniel Fernandez Perez
    volcano_plot <- ggplot(data = df, aes(x=log2FoldChange, y=-log10(padj), colour=DEG) ) +
                        geom_point(alpha=0.7, size=2) +

                        annotate("text", label = sum(df$DEG == "Upregulated"), color = "red", y = 0, x = 8, 
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
    print(volcano_plot)
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


# function to create the MA plots for genes and TEs

#read the table which tells you for each TE its family
if (snakemake@params[["organims"]] == "human"){
  families <- read.table("resources/TEs_per_family_hg38.txt", sep = "\t", col.names = c("TE","Family"))
} else if (snakemake@params[["organims"]] == "mouse"){
  families <- read.table("resources/TEs_per_family_mm10.txt", sep = "\t", col.names = c("TE","Family"))
}

#########################################################################
#--------------------------- MA plot  ----------------------------------#
#########################################################################

to_plot <- data.frame("Log2_mean_expression" = log2(table_tes$baseMean+1),
                      "Log2_FoldChange" =table_tes$log2FoldChange,
                      "Significance"=table_tes$DEG) 

rownames(to_plot) <- table_tes$Geneid

to_plot$Family <- families[match(rownames(to_plot),families$TE),"Family"] # match the two files

to_plot[which(is.na(to_plot$Family)),"Family"] <- "Other" #in case there are NAS set a label

#set family for each TE
families_counts <- table(to_plot$Family)
to_keep <- names(which(families_counts > 10))
to_plot$Family <- ifelse(to_plot$Family %in% to_keep, to_plot$Family, "Other")

to_plot <- to_plot %>% rownames_to_column(var = "Name")
#set threshold for minim expression
mean_thresh <- mean(to_plot$Log2_mean_expression)-sd(to_plot$Log2_mean_expression)

#plot
gplot_total <- ggplot(data = to_plot, aes(x=Log2_mean_expression, 
                                          y = Log2_FoldChange, 
                                          color = Family,
                                          label = Name))+
  geom_point(stat = "identity", color ="gray50")+
  geom_hline(yintercept=0)+
  geom_hline(aes(yintercept = log2FC),  color="red", linetype="dashed")+
  geom_hline(aes(yintercept = -log2FC),  color="red", linetype="dashed")+
  geom_vline(xintercept = mean_thresh, color ="red", linetype="dashed")+
  theme_classic()+
  geom_text_repel(aes(label=ifelse(abs(Log2_FoldChange)>log2FC & Log2_mean_expression > mean_thresh, as.character(Name),'')))+
  labs(title = paste0("TEs - ", plot_title))+
  xlab("Log2 (mean expression+1)")+
  ylab("Log2 (FC)")



pdf(snakemake@output[["ma_plot_TEs"]])
print(gplot_total)
dev.off()