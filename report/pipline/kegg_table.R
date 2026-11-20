args <- commandArgs()
group <- args[5]
#group <- 'WB_vs_SB'

path <- 'Report/6.Function/KEGG/'
#KEGG
file <- list.files(path,pattern = '.xls$')
file_name <- file[grep(paste0("^",group,'.kegg_enrichment.xls$'),file)]
list <- read.table(paste(path,file_name,sep=''),sep = '\t',quote = '',check.names = F,header = T,stringsAsFactors = F)

tmp <- list[,c(2,3,4,4,9,10)]
colnames(tmp)[4] <- 'pathview'
#tmp <- tmp[order(tmp[,1]),]

out <- tmp
out[,3] <- paste('<a href=\'',path,group,'_heatmap','/',gsub(' ','_',out[,3]),'.pdf','\' target="view_window">',out[,3],'</a>',sep = '')
out[,4] <- paste('<a href=\'',path,group,'_pathview','/',gsub(' ','_',out[,4]),'.pathview.png','\' target="view_window">',out[,4],'</a>',sep = '')
write.table(out,'kegg_table.txt',sep = '\t',quote = F,row.names = F,col.names = T)


