message("Loading libraries...\n")
suppressPackageStartupMessages(require(GenomicAlignments))
suppressPackageStartupMessages(require(rtracklayer))
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(GenomicRanges))
suppressPackageStartupMessages(require(IRanges))
suppressPackageStartupMessages(require(dplyr))



bamfiles <- snakemake@input[["bam_second_pass"]]

ifelse(length(bamfiles) == 0, stop("No bam files detected...\n"), no = 1)

#prepare gtf files
message("Reading the TEs gtf file...\n")
gtf_file_TE <- import(snakemake@params[["gtf_TEs"]])

#collapse all the regions from the same id
list_TE_is <- split(gtf_file_TE, gtf_file_TE$gene_id)

list_bams <- BamFileList(bamfiles)

#count reads on TEs
message("Count the reads mapped on the TEs...\n")

# inter.feat  = F: count multimapp
reads_counts_TE  <- summarizeOverlaps(list_TE_is, list_bams, mode = "IntersectionStrict", ignore.strand = F, inter.feature = F) 
tab_final_TE <- assay(reads_counts_TE)

colnames(tab_final_TE) <- sub(".bam", "", colnames(tab_final_TE))

tab_final_TE <- rownames_to_column(as.data.frame(tab_final_TE),var = "TE_id" )

write.table(tab_final_TE, file = snakemake@output[["tes_expression_counts"]], quote = F, sep = "\t", row.names = F )

message("Merge TEs and Gene expression...\n")

#merge GENE counts and TEs counts in one single table
if ( file.exists(snakemake@input[["gene_count_expression"]])){
    gene_counts <- read.table(snakemake@input[["gene_count_expression"]], header = T, sep = "\t", check.names = F)
    #check that the columns are exactly the same in both tables
    if ( all(colnames(gene_counts)[2:ncol(gene_counts)] == colnames(tab_final_TE)[2:ncol(tab_final_TE)]) ){
        colnames(gene_counts)[1] <- "Gene_TE_id"; colnames(tab_final_TE)[1] <- "Gene_TE_id"
        TEs_genes_table <- rbind(gene_counts, as.data.frame(tab_final_TE))
        write.table(TEs_genes_table,file = snakemake@output[["gene_tes_expression_counts"]], quote = F, sep = "\t", row.names = F)
    } else{
        message("Different column names between gene and TEs tables...\n")
    }
} else {
    message("The gene expression table does not exist...\n")
}