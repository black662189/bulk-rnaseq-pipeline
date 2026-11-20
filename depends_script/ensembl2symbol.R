count <- read.table("gene_count_matrix.csv",sep=",",header=T,check.names=F)
ensembl2symbol <- read.table("ensembl2symbol.txt",sep="\t",header=F,check.names=F)
colnames(ensembl2symbol) <- c("gene_id","gene")
merge <- merge(count,ensembl2symbol,by="gene_id")
merge2 <- merge[,c(ncol(merge),2:(ncol(merge)-1))]
write.csv(merge2,"gene_symbol_count_matrix.csv",row.names=F,quote=F)

fpkm <- read.table("all.ensembl.fpkm.xls",sep="\t",header=T,check.names=F)
merge <- merge(fpkm,ensembl2symbol,by="gene_id")
merge2 <- merge[,c(ncol(merge),2:(ncol(merge)-1))]
write.table(merge2,"all.genes.fpkm.xls",sep="\t",row.names=F,quote=F)

tpm <- read.table("all.ensembl.tpm.xls",sep="\t",header=T,check.names=F)
merge <- merge(fpkm,ensembl2symbol,by="gene_id")
merge2 <- merge[,c(ncol(merge),2:(ncol(merge)-1))]
write.table(merge2,"all.genes.tpm.xls",sep="\t",row.names=F,quote=F)
