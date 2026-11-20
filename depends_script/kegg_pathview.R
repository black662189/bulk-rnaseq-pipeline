args <- commandArgs()

library("pathview")

DEGs <- read.table(args[5],header=T,quote="",sep="\t",check.names=F)
gene_list <- as.vector(DEGs$log2FoldChange)
names(gene_list) <- DEGs$gene_id
ENS <- names(gene_list) %in% names(gene_list)[grep("ENS",names(gene_list))]
gene_list <- gene_list[!ENS]

kegg_list <- read.csv(args[6],header=T,sep="\t",check.names=F)
k <- nrow(kegg_list)
if(k==0){
q()
}
for(i in 1:k){
pid <- gsub("ko",args[7],kegg_list[i,]$pathway)
try(pathview(gene.data=gene_list,pathway.id=pid,species=args[7],limit=list(gene=3,cpd=1),gene.idtype="symbol"),silent=T)
name <- as.character(kegg_list[i,]$level3)
name1 <- gsub("/","",name)
name1 <- gsub(" ","_",name1)
name1 <- gsub("_-_","-",name1)
name1 <- gsub("[()]","",name1)
command1 <- paste("mv ",kegg_list[i,]$pathway,".pathview.png ",name1,".pathview.png",sep="")
command2 <- paste("mv ",kegg_list[i,]$pathway,".png ",name1,".png",sep="")
command3 <- paste("mv ",kegg_list[i,]$pathway,".xml ",name1,".xml",sep="")
system(command1)
system(command2)
system(command3)
}
