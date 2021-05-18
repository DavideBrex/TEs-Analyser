message("Loading libraries...\n")
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(IRanges))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(DESeq2))
suppressPackageStartupMessages(require(IHW))
suppressPackageStartupMessages(require(apeglm))
suppressPackageStartupMessages(require(EnhancedVolcano))
#set log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


#set the contrast
contrast <- paste("condition_", snakemake@params[["contrast"]][1], "_vs_", snakemake@params[["contrast"]][2], sep = "")

#load the deseq2 objects
for (i in 1:3){
  dds <- readRDS(snakemake@input[[i]])
  
  # Set as reference for the GLM model the sample that is going to act as a control in the contrast
  dds$condition <- relevel(dds$condition, snakemake@params[["contrast"]][2])
  dds <- nbinomWaldTest(dds)
  
  # Apply IHW to weight p-values based on baseMean: https://bioconductor.org/packages/release/bioc/vignettes/IHW/inst/doc/introduction_to_ihw.html
  res <- results(dds, name = contrast, filterFun = ihw, alpha = 0.05)
  
  #shrink fold changes for low expression genes
  res <- lfcShrink(dds = dds, coef = contrast, res = res, type = "apeglm")

  # Put nice the DEG table: sort by p-value, Geneid to column, round to 2 deciamls all columns except pvalues
  res <- res[order(res$padj),]
  pval            <- res$pvalue
  padjust         <- res$padj
  res.filt        <- as.data.frame(res) %>% rownames_to_column(var = "Geneid") %>% round_df(2)
  res.filt$pvalue <- pval
  res.filt$padj   <- padjust
  
  # Transform padj values that are 1 to NA (just to make nicer the volcano),
  # this was the default behaviour of DESeq2, but IHW doesn't do it
  res.filt <- res.filt %>% mutate(padj = ifelse(padj == 1, NA, padj))
  write.table( res.filt, file = snakemake@output[[i]], sep = "\t", quote = F, row.names = F ) 
}

