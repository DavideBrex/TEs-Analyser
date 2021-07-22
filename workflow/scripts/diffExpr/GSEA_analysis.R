#set log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(ReactomePA))
suppressPackageStartupMessages(library(org.Mm.eg.db))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(org.Hs.eg.db))
source("workflow/scripts/general_functions.R")


# Perform the Gene Set Enrichment Analysis

# Read table
table_genes <- read.delim(snakemake@input[[1]])


# Set mouse or human references for the databases
if (snakemake@params[["genome"]] == "mouse") { 
  kegg.genome <- "mmu"
  pa.genome   <- "mouse"
  db          <- "org.Mm.eg.db"
} else if (snakemake@params[["genome"]] == "human") {
    kegg.genome <- "hsa"
    pa.genome   <- "human"
    db          <- "org.Hs.eg.db"
}


# Get UP and DOWN-regulated 
UP <- table_genes %>% 
  dplyr::filter(DEG == "Upregulated") %>% 
  dplyr::select(Geneid) %>% 
  pull %>%
  { gsub("\\..*","", .)} %>%
  bitr(fromType = "SYMBOL",toType = c("ENTREZID"), OrgDb = db) #for ensembl set fromType = ENSEMBL

DWN <- table_genes %>%  
  dplyr::filter(DEG == "Downregulated") %>% 
  dplyr::select(Geneid) %>% 
  pull %>%
  { gsub("\\..*","",.)} %>%
  bitr(fromType = "SYMBOL",toType = c("ENTREZID"), OrgDb = db) #for ensembl set fromType = ENSEMBL


#perform analysis
UP.go      <- goEnrichment(UP, ont = "BP", db = db)
UP.kegg    <- KEGGenrichment(UP, org = kegg.genome, db = db)
UP.pa      <- PAenrichment(UP, org = pa.genome)
DWN.go     <- goEnrichment(DWN, ont = "BP", db = db)
DWN.kegg   <- KEGGenrichment(DWN, org = kegg.genome, db = db)
DWN.pa     <- PAenrichment(DWN, org = pa.genome)
GSEA.hall  <- GSEA_enrichment(table_genes, "resources/h.all.v6.2.symbols.gmt")
GSEA.c2all <- GSEA_enrichment(table_genes, "resources/c2.all.v6.2.symbols.gmt")
GSEA.c3tft <- GSEA_enrichment(table_genes, "resources/c3.tft.v6.2.symbols.gmt")


# ------ save outputs to xlsx ---------
list_of_datasets <- list( "GO Upregulated"        = UP.go@result %>% filter(p.adjust < 0.05) %>% add_row(), 
                          "GO Downregulated"      = DWN.go@result %>% filter(p.adjust < 0.05) %>% add_row(),
                          "KEGG Upregulated"      = UP.kegg@result %>% filter(p.adjust < 0.05) %>% add_row(),
                          "KEGG Downregulated"    = DWN.kegg@result %>% filter(p.adjust < 0.05) %>% add_row(),
                          "Reactome Upregulated"  = UP.pa@result %>% filter(p.adjust < 0.05) %>% add_row(),
                          "Reactome Downregulaed" = DWN.pa@result %>% filter(p.adjust < 0.05) %>% add_row(),
                          "GSEA Hallmarks"        = GSEA.hall %>% filter(padj < 0.25) %>% dplyr::select(-c(leadingEdge,nMoreExtreme)) %>% arrange(pval) %>% add_row(),
                          "GSEA c2all"            = GSEA.c2all %>% filter(padj < 0.25) %>% dplyr::select(-c(leadingEdge,nMoreExtreme)) %>% arrange(pval) %>% add_row(),
                          "GSEA c3tft"            = GSEA.c3tft %>% filter(padj < 0.25) %>% dplyr::select(-c(leadingEdge,nMoreExtreme)) %>% arrange(pval) %>% add_row())

write.xlsx(list_of_datasets, file = snakemake@output[["enrichments"]])
