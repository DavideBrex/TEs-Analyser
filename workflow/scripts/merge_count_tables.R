#set the log
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# Merge the read count tables for each sample in a single table #
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(rtracklayer))

# which coulmn to pick for each sample count matrix (see config.yaml for an explanation)
column_to_pick <- snakemake@params[["col_to_pick"]]

message("Reading the tables...\n")
#read all the files
datalist = lapply(snakemake@input, function(x){
    df <- read.table(file=x,header=F, sep = '\t')
    df <- df[,c(1,column_to_pick+1)]
    df <- df[grepl("ENS", df$V1),] #remove rows not starting with ENSEMBL identifier
    colnames(df) <- c("Gene_id", head(tail(strsplit(x, split = "/")[[1]], n=2), n=1))
    return(df)
})

message("Merging the tables...\n")  
#merge the tables
final_tab <- Reduce(function(x,y) {merge(x,y, by= 'Gene_id')}, datalist)

#now we convert from gene ensembl id to gene name (for those that have one)


gtf_file <- import(snakemake@params[["gtf_file"]]) #load gtf file
gtf_file <- as.data.frame(gtf_file)
gtf_file <- gtf_file[,c("gene_id","gene_name")]
gtf_file <- distinct(gtf_file) #select only one gene-id:gene_name tuple

#replace the gene name
final_tab$Gene_id <- gtf_file[match(final_tab$Gene_id, gtf_file$gene_id),"gene_name"]
message("Storing the table...\n")  
#store the gene count table
write.table( final_tab, file=snakemake@output[[1]], sep = "\t", quote = F, row.names = F )

message("All done!\n")  
