library(DESeq2)


args <- commandArgs()
data <- read.csv(args[5],header=T,sep=",",check.names=F)
#count matrix
countData <- data[,-1]
keep <- apply(countData,1,sum) > 0
countData <- countData[keep,]
keepData <- data[keep,]
#meta data 
condition <- factor(gsub('-[0-9]+','',colnames(data[,-1])))
dds <- DESeqDataSetFromMatrix(countData=countData,colData=DataFrame(condition),design=~condition)
dds2 <- DESeq(dds)
normalized_count <- counts(dds2,normalized = T)
normalized_count <- cbind(gene=as.character(keepData[,1]),normalized_count)
write.table(normalized_count,'all.genes.normalized.counts.xls',sep='\t',quote = F,row.names = F)
