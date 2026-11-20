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

RemoveGene <- function(exp,ratio=1.0){
  gua <- exp[,-1]
  keep <- c()
    for(i in 1:nrow(gua)){
			    zero_num <-sum(gua[i,]==0)
      if(zero_num>=ratio*ncol(gua)){
			        keep <- append(keep,i)
	      }
	    }
    out <-exp[-keep,]
	  return(out)
}

nodup_exp = myfunction(expr,"gene")

exp_result = RemoveGene(nodup_exp)

write.table(exp_result,output,sep="\t",row.names=F,quote=F)
