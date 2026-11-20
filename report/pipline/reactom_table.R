args <- commandArgs()
group <- args[5]
#group <- 'WB_vs_SB'

path <- 'Report/6.Function/Reactome/'
#reactom
file <- list.files(path,pattern = '.xls$')
file_name <- file[grep(paste0("^",group,'.reactome_enrichment.xls$'),file)]
list <- read.table(paste(path,file_name,sep=''),sep = '\t',quote = '',check.names = F,header = T,stringsAsFactors = F)

tmp <- list[,c(2,5,7,9)]
#tmp <- tmp[order(tmp[,1]),]
out <- tmp
out[,1] <- gsub('/','',out[,1])
out[,1] <- paste('<a href=\'',path,group,'_heatmap','/',gsub(' ','_',out[,1]),'.pdf','\' target="view_window">',out[,1],'</a>',sep = '')
write.table(out,'reactom_table.txt',sep = '\t',quote = F,row.names = F,col.names = T)


