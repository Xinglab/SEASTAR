#### FDR.R
#Using FDR.R to calculate the adjusted P-Value

argv <- commandArgs(TRUE)
path <- argv[1]


  input <- as.matrix(read.delim(paste(path,"/tmp/RMATS/rMATS_Result.txt",sep=""), header = TRUE))
  input <- input[order(input[,1]),]
  input[,"FDR"] <- p.adjust(as.numeric(input[,"PValue"]),method="fdr")
  ann<-read.table(file=paste(path,"/FilteredNrtss.annotation",sep=""), header=TRUE)
  ann<-ann[,1:10]
  colnames(ann)[1]<-'ID'
  input<-merge(ann,input,by='ID')
  input<-input[order(as.numeric(as.vector(input[,"FDR"]))),]
  write.table(input, paste(path,"/tmp/RMATS/rMATS_Result.txt",sep = ""),sep="\t", row.names=FALSE, col.name=TRUE, quote = FALSE)
  write.table(input, paste(path,"/rMATS_Result.txt",sep = ""),sep="\t", row.names=FALSE, col.name=TRUE, quote = FALSE)

#rMATS_Result.txt: Annotation of first exons and their testing results (a standard output format as the rMATS documentation) with FDR
