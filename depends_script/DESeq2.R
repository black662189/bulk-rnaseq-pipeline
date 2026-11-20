library(DESeq2)


args <- commandArgs()
data <- read.csv(args[5],header=T,sep=",",check.names=F)
compare <- read.table(args[6],header=F,sep="\t",quote = '',check.names=F,stringsAsFactors = F)
foldchange <- as.numeric(args[7])
qvalue <- as.numeric(args[8])

data$sum <- apply(data[,-1],1,sum)
data <- data[rev(order(data$sum)),]
data <- data[!duplicated(data[,1]),]
data <- data[,-ncol(data)]

#print (foldchange)
#count matrix
countData <- data[,-1]
#keep <- apply(countData,1,sum) >= 10
### modify by yueyao at 20191211
keep <- apply(countData,1,sum) >0

countData <- countData[keep,]
keepData <- data[keep,]

#meta data 
condition <- factor(gsub('-[0-9]+','',colnames(data[,-1])))
dds <- DESeqDataSetFromMatrix(countData=countData,colData=DataFrame(condition),design=~condition)
dds2 <- DESeq(dds)
normalized_count <- counts(dds2,normalized = T)
normalized_count <- cbind(gene=as.character(keepData[,1]),normalized_count)
write.table(normalized_count,'all.genes.normalized.counts.xls',sep='\t',quote = F,row.names = F)
#result
for(c in 1:nrow(compare)){
  ctrl <- compare[c,1]
  treat <- compare[c,3]
  ctrl_columns <- grep(paste0("^",ctrl,"-"),colnames(keepData))
  treat_columns <- grep(paste0("^",treat,"-"),colnames(keepData))
  #raw counts
  gua <- keepData[,c(ctrl_columns,treat_columns)]
  #normalized counts
  gua1 <- normalized_count[,c(ctrl_columns,treat_columns)]
  colnames(gua1) <- paste0(colnames(gua1),'(norm)')
  res <- results(dds2,contrast=c("condition",treat,ctrl))
  res2 <- cbind(gene_id=keepData[,1],total_count=apply(gua,1,sum),gua,gua1,as.matrix(res))
  
  resOrdered <- res2[order(res2$padj),]
  write.table(resOrdered,file=paste(ctrl,"_vs_",treat,".total_genes.xls",sep=""),sep="\t",quote=F,row.names=F)
  
  DEGs_FC <- abs(resOrdered$log2FoldChange) >= log2(foldchange) & resOrdered$padj < qvalue & !is.na(resOrdered$padj)
  write.table(resOrdered[DEGs_FC,],file=paste(ctrl,"_vs_",treat,".differentially_expressed_genes.xls",sep=""),sep="\t",quote=F,row.names=F)
}
