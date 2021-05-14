# Merge the read count tables for each sample in a single table #

# which coulmn to pick for each sample count matrix (see config.yaml for an explanation)
column_to_pick <- snakemake@params[["col_to_pick"]]

#read all the files
datalist = lapply(snakemake@input, function(x){
    df <- read.table(file=x,header=F, sep = '\t')
    df <- df[,c(1,column_to_pick+1)]
    df <- df[grepl("ENS", df$V1),] #remove rows not starting with ENSEMBL identifier
    colnames(df) <- c("Gene_id", head(tail(strsplit(x, split = "/")[[1]], n=2), n=1))
    return(df)
})

#merge the tables
final_tab <- Reduce(function(x,y) {merge(x,y, by= 'Gene_id')}, datalist)

#store the result
write.table( final_tab, file=snakemake@output[[1]], sep = "\t", quote = F, row.names = F )