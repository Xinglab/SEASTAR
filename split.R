#### split.R
#Using split.R to segment the 5'UTR into two regions based on change-point model

argv <- commandArgs(TRUE)
inpath <- argv[1]
path <- argv[2]

#inpath <- "E:/AFE/iPS/utr_result_4,5,6_1,2,3_All_Prediction_Results.txt"
#path <- "E:/AFE/iPS/"


  input <- read.delim(paste(inpath,sep=""), header = TRUE)
  input <- input[order(input[,"Gene"]),]
  ann <- read.table(file=paste(path,"/FilteredNrtss.annotation",sep=""), header=TRUE)
  ann <- merge(ann,input[,1:3],by.x='tss_id', by.y="Gene")
  left <- ann[,c("tss_chr", "tss_ss", "Predicted_Proximal_5UTR", "tss_id", "in.out", "strand")]
  left[,"tss_id"] <- paste(left[,"tss_id"],".1",sep="")
  right <- ann[,c("tss_chr", "Predicted_Proximal_5UTR", "tss_es", "tss_id", "in.out", "strand")]
  right[,"Predicted_Proximal_5UTR"] <- right[,"Predicted_Proximal_5UTR"] + 1
  right[,"tss_id"] <- paste(right[,"tss_id"],".2",sep="")
  colnames(left) <- ""
  colnames(right) <- ""
  output <- rbind(left, right)
  output <- output[order(output[,4]),]
  write.table(ann, paste(path,"/utr.annotation",sep = ""),sep="\t", row.names=FALSE, col.name=TRUE, quote = FALSE)
  write.table(output, paste(path,"/utr_split.bed",sep = ""),sep="\t", row.names=FALSE, col.name=F, quote = FALSE)

#utr.annotation: Annotation file of utr adding with fit value and changed point predicted by Dapars prediction model
#utr_split.bed: bed region of utr split into two regions
