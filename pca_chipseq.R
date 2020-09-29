library(ggplot2)
library(ggrepel)
library(stringr)
library(DESeq2)
library(edgeR)
coverage <- read.csv('/media/shaoqizhu/easystore/hdac/bampe/coverage.bed',sep='\t',header = F)
name <- read.csv('/media/shaoqizhu/easystore/hdac/bampe/name.txt',sep='\t',header = F)
merge <- read.csv('/media/shaoqizhu/easystore/hdac/bampe/merged.bed',sep='\t',header = F)
names(merge)<-c('chrome','start','end','idx')

RPK <- coverage[,2:11]/coverage[,1]*10^3
RPKM <- as.data.frame(t(t(RPK)/apply(RPK, 2, sum)*10^7))

pca <- as.data.frame(prcomp(t(coverage[,2:11]),scale.=T)$x)
pca <- as.data.frame(prcomp(t(RPKM[,]),scale.=T)$x)
pca$name <- name$V1
pca$exp <- paste(str_split_fixed(name$V1,'_',3)[,1],str_split_fixed(name$V1,'_',3)[,2],sep = '_')
pca$exp[5:6] <- "na"

ggplot(data = pca,aes(x=pca$PC1, y=pca$PC2, label=name, colour = exp)) + 
  geom_point() + geom_text_repel(aes(label=name), size = 4, hjust=0, vjust=-1, angle=0)

PCA <- prcomp(t(coverage[2:11]),scale.=T)
PCA <- prcomp(t(RPKM[,]),scale.=T)
(PCA$sdev/sum(PCA$sdev))[1:2]

diff <- numeric(); N = length(unique(pca$exp));Name <- character()
for (i in 1:(N-1)){
  for (j in (i+1):N){
    condition1 = which(pca$exp==unique(pca$exp)[i])
    condition2 = which(pca$exp==unique(pca$exp)[j])
    coldata <- as.data.frame(pca$exp[c(condition1,condition2)])
    colnames(coldata) = 'condition'
    rownames(coldata) = as.character(names(coverage)[c(condition1,condition2)+1])
    dds <- DESeqDataSetFromMatrix(countData = coverage[,c(condition1,condition2)+1],
                              colData = coldata,
                              design = ~ condition)
    dds <- DESeq(dds,fitType = 'local')
    res <- results(dds)
    print(c(condition1,condition2)+1)
    #diff_idx <- which(res$padj < 0.01 & abs(res$log2FoldChange) > log2(3))
    #diff <- c(diff,diff_idx)
    diff = cbind(diff,res$log2FoldChange,res$padj)
    Name <- c(Name,paste(unique(pca$exp)[j],'/',unique(pca$exp)[i],'(log2FC)',sep=''),paste(unique(pca$exp)[j],'/',unique(pca$exp)[i],'(padj)',sep=''))
    #diff <- c(diff,diff_0)
    #change <- cbind(merge[diff_idx,],res$log2FoldChange[diff_idx])
    #write.table(change,paste('/media/shaoqizhu/easystore/hdac/bampe/diff_peaks/',unique(pca$exp)[i],'-',unique(pca$exp)[j],'.bed',sep=''),
    #            quote = F,sep = '\t',row.names = F,col.names = F)
  }
}

diff <- as.data.frame(diff)
names(diff) <- Name; names(merge) <- c('chrome','start','end','idx');names(RPKM) <- name$V1
write.table(cbind(merge,round(RPKM,4),round(diff[,c(3,4,5,6,15,16)],4)),paste('/media/shaoqizhu/easystore/hdac/bampe/','merge_peak_diff.bed',sep=''),
            quote = F,sep = '\t',row.names = F,col.names = T)

write.table(merge[unique(diff),],paste('/media/shaoqizhu/easystore/hdac/pca_no_dkmyc/','merge_diff_peak','.bed',sep=''),
                        quote = F,sep = '\t',row.names = F,col.names = F)

diff <- numeric(); N = length(unique(pca$exp))
for (i in 1:(N-1)){
  for (j in (i+1):N){
group <- factor(c(1,1,2,2))
y = DGEList(counts = coverage[,c(2*i,2*i+1,2*j,2*j+1)],group = group)
print(c(2*i,2*i+1,2*j,2*j+1))
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y<-estimateTagwiseDisp(y)
#design <- model.matrix(~group)
#fit <- glmQLFit(y, design)
#qlf <- glmQLFTest(fit, coef=2)
et <- exactTest(y)
diff=rbind(diff,length(which(et$table$PValue < 0.01 & abs(et$table$logFC) > log2(3))))
  }
}


A = read.csv('/media/shaoqizhu/easystore/hdac/pca_all/WT_D0-WT_D1.bed',sep='\t',header = F)
B = read.csv('/media/shaoqizhu/easystore/hdac/pca_all/dKO_D1-WT_D1.bed',sep='\t',header = F)
write.table(merge[intersect(A$V4[A$V5>0],B$V4[B$V5<0]),],'/media/shaoqizhu/easystore/hdac/pca_all/WT_D1_unique.bed',
            quote = F,sep = '\t',row.names = F,col.names = F)

dKO_D0_1_length <- read_table2('/media/shaoqizhu/easystore/hdac/bampe/bedpe/dKO_D1_2_length.bed',col_names = F)
hist(dKO_D0_1_length$X1,40,main)
max(dKO_D0_1_length$X1)
