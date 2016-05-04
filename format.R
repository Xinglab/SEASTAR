#### format.R
#Using format.R to estimate the the coverage and effective length of each first exon in one AFE event preparing for the estimation of PSI

argv <- commandArgs(TRUE)
path <- argv[1]
length1 <- argv[2]


#######################  identify the differential ratio of tss by multi test final ##############################################
#### run by rMATs

HBM <- read.delim(paste(path,"/tmp/table/table.cov",sep=''), header = FALSE)
##find gene_id coresponding tss_id (non-redundant tss), remove the CS and PTL, only remain the HBM and NSVD
table <- HBM[,-1]
rownames(table) <- HBM[,1]
colnames(table) <- 1:ncol(table)
table <- table[order(as.vector(HBM[,1])),]

tss <- read.delim(paste(path,"/tmp/tsgtf/nrtss_filter.annotation",sep=''), header=TRUE)
nrtss <- as.matrix(read.table(paste(path,"/tmp/tsgtf/nrtss.annotation",sep=''),header = TRUE, sep = "\t",quote = ""))
nrtss <- cbind(nrtss,tss$result)
colnames(nrtss)[ncol(nrtss)] <- "result"

table <- table[order(nrtss[,3],nrtss[,4],nrtss[,1]),] ##multi tss, add, Cuff gene
nrtss <- nrtss[order(nrtss[,3],nrtss[,4],nrtss[,1]),] ##multi tss, add, Cuff gene
table <- table[nrtss[,"result"]=="TRUE",] #filter
nrtss <- nrtss[nrtss[,"result"]=="TRUE",] #filter

##rMATS input examples
##calculate the length for each and the total count for each sample and each gene
#ID  		IC_SAMPLE_1	SC_SAMPLE_1	IC_SAMPLE_2	SC_SAMPLE_2	IncFormLen	SkipFormLen
#chr3:9422780	1,0		6,4		0,0,1		1,1,1		100		100
#chr3:9426702	0,0		9,15		1,0,0		18,14,21	100		100
table <- table*100/(as.numeric(nrtss[,9])-as.numeric(nrtss[,8])+1) #normalize by tss length
table <- table[((as.numeric(nrtss[,9])-as.numeric(nrtss[,8]))>=3),] #not short than 3
nrtss <- nrtss[((as.numeric(nrtss[,9])-as.numeric(nrtss[,8]))>=3),] #not short than 3

total <- table
for (i in 2:nrow(nrtss)) {
  if (nrtss[i-1,3]==nrtss[i,3]) { ##multi tss, add, Cuff gene
    total[i,] <- total[i-1,] + total[i,]
  }
}
for (i in nrow(nrtss):2) {
  if (nrtss[i-1,3]==nrtss[i,3]) { ##multi tss, add, Cuff gene
    total[i-1,] <- total[i,]
  }
}
tmp <- matrix(0,nrow(total),7)
colnames(tmp) <- c("ID","IC_SAMPLE_1","SC_SAMPLE_1","IC_SAMPLE_2","SC_SAMPLE_2","IncFormLen","SkipFormLen")
tmp[,1] <- nrtss[,1]
tmp[,6] <- as.numeric(nrtss[,9])-as.numeric(nrtss[,8])+1 ## use tss length
tmp[,7] <- as.numeric(nrtss[,9])-as.numeric(nrtss[,8])+1 ## use tss length

table <- table*(as.numeric(nrtss[,9])-as.numeric(nrtss[,8])+1)/100 #recover by tss length
total <- total*(as.numeric(nrtss[,9])-as.numeric(nrtss[,8])+1)/100 #recover by tss length

table <- round(table)
total <- round(total)
rest <- total-table
rest <- round(rest)
    write.table(rest, paste(path,"/tmp/table/rest_multitest.cov",sep=""),sep="\t", row.names=TRUE, col.name=TRUE, quote = FALSE)
    write.table(total, paste(path,"/tmp/table/total_multitest.cov",sep=""),sep="\t", row.names=TRUE, col.name=TRUE, quote = FALSE)
    write.table(table, paste(path,"/tmp/table/table_multitest.cov",sep=""),sep="\t", row.names=TRUE, col.name=TRUE, quote = FALSE)


for (i in 1:nrow(tmp)) {
  tmp[i,2] <- paste(table[i,c(1:as.numeric(length1))],collapse=",")
  tmp[i,3] <- paste(rest[i,c(1:as.numeric(length1))],collapse=",")
  tmp[i,4] <- paste(table[i,c((as.numeric(length1)+1):ncol(table))],collapse=",")
  tmp[i,5] <- paste(rest[i,c((as.numeric(length1)+1):ncol(table))],collapse=",")
}
#Checking the input to satisfy requirements of rMATS
output <- tmp[((rowSums(as.matrix(total[,c(1:as.numeric(length1))],nrow(total)))!=0)+(rowSums(as.matrix(total[,c((as.numeric(length1)+1):ncol(table))],nrow(total)))!=0)==2),]  ##not divide by 0 calculating PSI
write.table(output, paste(path,"/tmp/table/RMATS_AB_input.txt",sep=""),sep="\t", row.names=FALSE, col.name=TRUE, quote = FALSE)

#Seprating the input into several smaller files for calculating in rMATS model
#To prevent the model hangs on or crashes with a large file as the input
output <- data.frame(output)
bin <- min(1000,nrow(output))  ##split into 1000 pieces or the min rows
for (j in 1:bin) {
  write.table(output[c((floor((nrow(output)*(j-1)/bin))+1):floor(nrow(output)*j/bin)),], paste(path,"/tmp/RMATS/RMATS_",j,".txt",sep = ""),sep="\t", row.names=FALSE, col.name=TRUE, quote = FALSE)
}

#RMATS_n.txt: the standard file format as the input of rMATS model