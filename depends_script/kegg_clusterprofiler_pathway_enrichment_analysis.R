#@Time:  15:20
#@Author:ldyyyy
#@File:kegg_clusterprofiler_pathway_enrichment_analysis.py
#@Software:PyCharm
#@Email:ldyyyy21@foxmail.com

args <- commandArgs()
degfile <- args[6]
species <- args[7]
groupname <- args[8]
del_human = T
gmtout = F
getko = F
outkopath = F
p = 0.05
fdr = 1
zs = 0
library(ggplot2)
library(RColorBrewer)
library(Orgdb)
library(clusterProfiler)
### 1.通路富集
deg <- read.table(degfile,header=T,sep="\t",quote="",stringsAsFactors=F)
gene <- deg$gene_id
entrezID2symbol <- AnnotationDbi::select(Orgdb,keys=as.character(gene),columns='ENTREZID',keytype='SYMBOL')
kegg  <- enrichKEGG(entrezID2symbol$ENTREZID,organism=species,pvalueCutoff=1,qvalueCutoff=1,keyType='kegg')
kegg <- setReadable(kegg, OrgDb = Orgdb, keyType="ENTREZID")
kegg_out <- data.frame(kegg)
kegg_out <- kegg_out[order(kegg_out$pvalue),]
kegg_out_filter <- kegg_out[kegg_out$pvalue<p,]
kegg_out_filter <- kegg_out_filter[kegg_out_filter$`p.adjust`<fdr,]
kegg_out_filter <- kegg_out_filter[kegg_out_filter$Count>=3,]
write.table(kegg_out,paste("unfiltered.",groupname,".kegg_enrichment.xls",sep=""),sep="\t",row.names=F,quote=F)
write.table(kegg_out_filter,paste(groupname,".kegg_enrichment.xls",sep=""),sep="\t",row.names=F,quote=F)

### 2.绘图
## 数据准备
pathway <- kegg_out_filter
if(nrow(pathway) ==0){
  cat(c("Warning ",groupname,".kegg_enrichment.xls is blank\n"))
  q()
}

if(nrow(pathway)>15) {
data <- pathway[1:15,]
} else {
data <- pathway
}



if(nrow(data)>10) {
  h <- 8*nrow(data)/60
} else if(nrow(data)<=5){
  h <- 1.2
}else{
  h <- 1.6
}

# shorten_names <- function(x, n_word=4, n_char=40){
#   shortname = x["Description"]
#   if (length(strsplit(shortname, " ")[[1]]) > n_word || (nchar(shortname) > n_char))
#   {
#     if (nchar(shortname) > n_char) shortname <- substr(shortname, 1, n_char)
#     y <- paste(paste(strsplit(shortname, " ")[[1]][1:min(length(strsplit(shortname," ")[[1]]), n_word)],
#                      collapse=" "), "...", sep="")
#     return(y)
#   }
#   else
#   {
#     return (shortname)
#   }
# }
#
# data$pretty_varname = apply(data,1,shorten_names)
max1<- max(nchar(c(as.character(data$Description))))
w<- max1*0.09+1.5


## 柱状图

bar <- ggplot(data,aes(y=Count,x=reorder(Description,-log10(pvalue)))) +
  geom_bar(aes(fill=-log10(pvalue)),stat="identity",width=0.8,position = position_dodge(0.9)) +
  scale_fill_gradientn(colours=c(rev(brewer.pal(11, "RdYlBu")[2:10]))) +
  #scale_fill_gradientn(colours=c(rev(brewer.pal(11, "Spectral")[1:5]))) +
  coord_flip() +labs(color=expression(-log10(pvalue)),y="Gene Number",x="Pathway") +
  theme_bw() +
  theme(
		title =element_text(size=6,face="plain",colour="black"),
		axis.title = element_text(size=6,face="plain",colour="black"),
		panel.grid=element_blank(),

		axis.line = element_line(size=0.25,colour = "black"),
		axis.text = element_text(color = "black",size = 6,face="plain"),

		legend.key.size=unit(0.3,'cm'),
		legend.title=element_text(size=6,colour = "black",face="plain"),
		legend.text = element_text(colour = 'black',face = 'plain',size=6)

)
ggsave(paste(groupname,".kegg_enrichment_bar.png",sep=""),bar,width=w,height=h)
ggsave(paste(groupname,".kegg_enrichment_bar.pdf",sep=""),bar,width=w,height=h)
#ggsave(paste(groupname,".kegg_enrichment_bar.svg",sep=""),bar,width=w,height=h)
svg_file=paste(groupname,".kegg_enrichment_bar.svg",sep="")
svg(svg_file,width=w,height=h)
plot(bar)
dev.off()


