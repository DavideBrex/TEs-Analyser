quantify_TEs_expr <- function(path_to_main, TEs_gtf) {
  #check that the files exist and are in the correct folder
  setwd(paste0(path_to_main,"AllBAMs/"))
  bamfiles <- list.files(pattern = "*tagged_reads.filtered*", recursive = T)
  ifelse(length(bamfiles) == 0, stop("No bam files detected...\n"), no = 1)
  message("Loading libraries...\n")
  suppressPackageStartupMessages(require(GenomicAlignments))
  suppressPackageStartupMessages(require(rtracklayer))
  suppressPackageStartupMessages(require(tidyverse))
  suppressPackageStartupMessages(require(GenomicRanges))
  suppressPackageStartupMessages(require(IRanges))
  suppressPackageStartupMessages(require(dplyr))
  #prepare gtf files
  message("Reading the TEs gtf file...\n")
  gtf_file_TE <- import(TEs_gtf)
  list_TE_is <- split(gtf_file_TE, gtf_file_TE$gene_id)
  
  list_bams <- BamFileList(bamfiles)
  print(bamfiles)
  #count reads on TEs
  message("Count the reads mapped on the TEs...\n")
  reads_counts_TE  <- summarizeOverlaps(list_TE_is, list_bams, mode = "IntersectionStrict",ignore.strand = F,
                                        inter.feature = F) # inter.feat  = F: count multimapp
  tab_final_TE <- assay(reads_counts_TE)
  colnames(tab_final_TE) <- sub(".tagged_reads.filtered.bam", "", colnames(tab_final_TE))
  tab_final_TE <- rownames_to_column(as.data.frame(tab_final_TE),var = "TE_id" )
  setwd(path_to_main)
  write.table(tab_final_TE, "TEs_expression_counts.txt", quote = F, sep = "\t", row.names = F )
  message("Merge TEs and Gene expression...\n")
  #merge GENE counts and TEs counts in one single table
  if ( file.exists("Reads_count_HTseq.txt")){
    gene_counts <- read.table("Reads_count_HTseq.txt", header = T, sep = "\t", check.names = F)
    #check that the columns are exactly the same in both tables
    if ( all(colnames(gene_counts)[2:ncol(gene_counts)] == colnames(tab_final_TE)[2:ncol(tab_final_TE)]) ){
      colnames(gene_counts)[1] <- "Gene_TE_id"; colnames(tab_final_TE)[1] <- "Gene_TE_id"
      TEs_genes_table <- rbind(gene_counts, as.data.frame(tab_final_TE))
      write.table(TEs_genes_table,"Reads_count_TEs_genes.txt", quote = F, sep = "\t", row.names = F)
    }
    else{
      message("Different column names between gene and TEs tables...\n")
    }
  } else {
    message("The gene expression table does not exist...\n")
  }
}