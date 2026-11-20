args <- commandArgs()
group <- args[5]
#group <- 'SB_vs_SP'

path <- 'Report/6.Function/GO/'
#biological_process
#count_png <- list.files(path,pattern = paste0(group,'.*.png$'))
#tmp <- gsub('.png','',count_png)
tmp <- list.files(path,pattern = paste0(group,'\\..*.xls$'))
for (i in tmp) {
  #提取过程名pb/cc/mf
  part <- gsub(paste(group,'_enrichment.xls',sep = '|'),'',i)
  part <- gsub('\\.','',part)

  list <- read.table(paste0(path,i),sep = '\t',quote = '',check.names = F,header = T,stringsAsFactors = F)
  tmp1 <- list[,c(2,5,7,9)]
  #tmp1 <- tmp1[order(tmp1[,1]),]
  out <- tmp1
  out$Description <- gsub('/','',out$Description)
  if (dim(out)[1]==0){
    tmp2 <- gsub('_enrichment','',unlist(strsplit(i,'\\.'))[2])
    write.table(out,paste0('GO_table_',tmp2,'.txt'),sep = '\t',quote = F,row.names = F,col.names = T)
  } else{
    out[,1] <- paste('<a href=\'',path,group,'_heatmap','/',part,'/',gsub(' ','_',out[,1]),'.pdf','\' target="view_window">',out[,1],'</a>',sep = '')
	tmp2 <- gsub('_enrichment','',unlist(strsplit(i,'\\.'))[2])
    write.table(out,paste0('GO_table_',tmp2,'.txt'),sep = '\t',quote = F,row.names = F,col.names = T)
  }
  
}
