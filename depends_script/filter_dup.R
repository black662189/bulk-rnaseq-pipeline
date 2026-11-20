args <- commandArgs()

expr <- read.table(args[5],sep="\t",header=T,check.names=F,quote="")
output <- args[6]

myfunction <- function(rt,IdColName){
  rt1 <- rt
  idColNum <- grep(IdColName,colnames(rt))
    mat <- as.matrix(rt[,-idColNum])
    rownames(mat) <- rt[,idColNum]
	  rt1$mean <- apply(mat, 1, sum)
	  aa <- rt1[order(rt1[,IdColName],rt1$mean,decreasing = T),]
	    aa <- aa[!duplicated(aa[,IdColName]),]
	    aa <- subset(aa,select=-c(mean))
		  return(aa)
}


nodup_exp = myfunction(expr,"gene")

write.table(nodup_exp,output,sep="\t",row.names=F,quote=F)
