exonLen <- read.table('/media/shaoqizhu/easystore/Annotation/mm10_exonLen.txt')
counts <- read.table('/media/shaoqizhu/easystore/qshan/counts.txt',row.names = 1)
name <- read.table('/media/shaoqizhu/easystore/qshan/name.txt')
names(counts) <- name$V1
counts <- counts[order(rownames(counts)),]
counts <- counts[-c(1:5),]
RPK <- counts/exonLen$V2*10^3
TPM <- as.data.frame(t(t(RPK)/apply(RPK, 2, sum)*10^6))
write.csv(TPM,'/media/shaoqizhu/easystore/qshan/TPM.csv')
    