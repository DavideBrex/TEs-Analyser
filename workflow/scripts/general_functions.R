
#Function to round all numeric colums of dataframe
round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}

########## ENRICHMENTS FUNCTIONS #################

goEnrichment <- function(df, ont = "BP", db = org.Mm.eg.db) {
  require(clusterProfiler)
  require(db)
  ego <- enrichGO(gene          = df$ENTREZID,
                  OrgDb         = db,
                  keyType       = 'ENTREZID',
                  ont           = ont,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.1,
                  readable = TRUE)
  message("GO")
  print(is.null(ego))
  return(ego)
}

KEGGenrichment <- function(df, org = "mmu", db = org.Mm.eg.db) {
  require(clusterProfiler)
  require(db)
  ekgg <- enrichKEGG(gene          = df$ENTREZID,
                     organism      = org,
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.1)
  ekgg <- DOSE::setReadable(ekgg, OrgDb = db, keyType="ENTREZID")
  message("KEGG")
  print(is.null(ekgg))
  return(ekgg)
}


PAenrichment <- function(df, org = "mouse") {
  require(clusterProfiler)
  require(ReactomePA)
  ePA <- enrichPathway(gene         = df$ENTREZID,
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.1,
                       organism     = org,
                       readable     = TRUE)
  message("Reactome")
  print(is.null(ePA))
  return(ePA)
}

GSEA_enrichment <- function(df, pathways.gmt) {
  require(clusterProfiler)
  pathways <- fgsea::gmtPathways(pathways.gmt)
  
  # Create a list containing a named vector (with genenames) of log2fc of each PcgfsKO vs WT
  geneList <- df$log2FoldChange
  names(geneList) <- toupper(df$Geneid)
  
  # Run GSEA algorithm
  fgseaRes <- fgsea::fgsea(pathways = pathways, 
                           stats   = geneList,
                           minSize = 15,
                           maxSize = 2000,
                           nperm   = 10000)
  return(fgseaRes)
}