args <- commandArgs()
#path <- args[5]
group <- args[5]
path <- 'Report/3.DEG/'
#group <- 'WB_vs_SB'
file <- list.files(path,pattern = '.differentially_expressed_genes.xls')
file <- file[grep(paste0('^',group,'.differentially_expressed_genes.xls$'),file)]
list <- read.table(paste(path,file,sep = ''),sep = '\t',quote = '',check.names = F,header = T)
#num <- dim(list)[2]
tmp <- list[,c('gene_id','log2FoldChange','pvalue','padj')]
out <- tmp[order(tmp[,2]),]
#out <- out[grep('[a-zA-Z]', out[,4],invert = T),]
down <- head(out,n=50)
up <- tail(out,n=50)
up <- up[order(up[,2],decreasing = T),]
write.table(up,'up_DEG_table.txt',sep = '\t',quote = F,row.names = F,col.names = T)
write.table(down,'down_DEG_table.txt',sep = '\t',quote = F,row.names = F,col.names = T)
