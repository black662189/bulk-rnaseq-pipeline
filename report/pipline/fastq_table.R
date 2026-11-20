args <- commandArgs()
#path <- args[5]
#file <- list.files(path)
file <- list.dirs('Report//1.QC//FastQC/')[-1]
#dir <- paste('Report//QC//FastQC//',file,'//',sep = '')
name <- gsub("Report//1.QC//FastQC//","",file)
tmp1 <- c()
tmp2 <- c()
file1 <- c()
file2 <- c()
for (i in file) {
  fastq <- list.files(i,pattern = '.html')
  tmp1 <- append(tmp1,fastq[1])
  tmp2 <- append(tmp2,fastq[2])
  file1 <- append(file1,sub('.html','',fastq[1]))
  file2 <- append(file2,sub('.html','',fastq[2]))  
}
tmp <- data.frame(file,tmp1,tmp2)
out <- data.frame(name,file1,file2)
out1 <- out

out1[,2] <- paste('<a href=\'',tmp[,1],'//',tmp[,2],'\' target="view_window">',out[,2],'</a>',sep = '')
out1[,3] <- paste('<a href=\'',tmp[,1],'//',tmp[,3],'\' target="view_window">',out[,3],'</a>',sep = '')

write.table(out1,'QC_table.txt',quote = F,sep = '\t',row.names = F,col.names = F)



