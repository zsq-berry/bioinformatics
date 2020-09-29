library(ggplot2)
library(ggrepel)
library(stringr)
library(edgeR)
library(grid)
library(factoextra)
library(cluster)

coverage <- read.csv('/media/shaoqizhu/easystore/Stat5/bam/coverage.bed',sep='\t',header = F)
merge <- read.csv('/media/shaoqizhu/easystore/Stat5/narrowPeak/merged.bed',sep='\t',header = F)
name <- read.csv('/media/shaoqizhu/easystore/Stat5/bam/name.txt',sep='\t',header = F)
names(coverage)[2:5] <- as.character(name$V1[1:4])
names(merge)<-c('chr','start','end','idx')

RPK <- coverage[,2:5]/coverage[,1]*10^3
TPM <- as.data.frame(t(t(RPK)/apply(RPK, 2, sum)*10^6))
coverage <- coverage[,-1]

pca <- as.data.frame(prcomp(t(TPM),scale.=T)$x)
pca$name <- name$V1[1:4];pca$sample <- name$V1[5:8]

ggplot(data = pca,aes(x=PC1, y=PC2, label=name, colour = sample)) + 
  geom_point() + geom_text_repel(aes(label=name), size = 4, hjust=0, vjust=-1, angle=0)


PCA <- prcomp(t(coverage[,]),scale.=T)
PCA <- prcomp(t(RPKM[,]),scale.=T)
(PCA$sdev/sum(PCA$sdev))[1:2]

diff <- numeric(); N = length(unique(pca$sample))
for (i in 1:(N-1)){
  for (j in (i+1):N){
    condition1 = which(pca$sample==unique(pca$sample)[i])
    condition2 = which(pca$sample==unique(pca$sample)[j])
    group <- factor(c(rep(i,length(condition1)),rep(j,length(condition2))))
    y = DGEList(counts = coverage[,c(condition1,condition2)],group = group)
    print(names(coverage[,c(condition1,condition2)]))
    y <- calcNormFactors(y)
    y <- estimateCommonDisp(y)
    y <- estimateTagwiseDisp(y)
    fit <- glmQLFit(y)
    qlf <- glmQLFTest(fit)
    fdr <- p.adjust(p = qlf$table$PValue,method = "BH",n = length(qlf$table$PValue))
    diff_idx <- which(fdr < 0.05 & abs(qlf$table$logFC) > log2(2))
    diff <- c(diff,diff_idx)
    change <- cbind(merge[diff_idx,],logFC=qlf$table$logFC[diff_idx],pvalue=qlf$table$PValue[diff_idx],fdr=fdr[diff_idx])
    #write.table(change,paste('/media/shaoqizhu/easystore/CD8-HP/diff_peaks/',unique(pca$exp)[i],'-',unique(pca$exp)[j],'.bed',sep=''),
    #           quote = F,sep = '\t',row.names = F,col.names = F)
  }
}
length(unique(diff))

stat5_diff <- change

hits <- findOverlaps(GRanges(stat5_diff),GRanges(diff_peak_kmeans))

stat5_DNase <- cbind(stat5_diff[hits@from,],diff_peak_kmeans[hits@to,c(4,5)])

sum(stat5_diff$logFC>0)

hits <- findOverlaps(GRanges(stat5_diff),GRanges(tcf1))

nrow(subset(stat5_diff[hits@from,],logFC<0))

#Stat5 dKO<WT
table(stat5_DNase$cluster[stat5_DNase$logFC>0])
#Stat5 dKO>WT
table(stat5_DNase$cluster[stat5_DNase$logFC<0])

hits <- findOverlaps(GRanges(stat5_diff),GRanges(mm9_deg_50kb))
stat5_deg <- cbind(stat5_diff[hits@from,],mm9_deg[hits@to,])
dist_tss <- sapply(1:nrow(stat5_deg), function(x){stat5_deg[x,9]-stat5_deg[x,2]})
stat5_deg$dist_tss <- dist_tss

# logFC<0: dko>wt; logFC>0: wt>dko
cbind(table(stat5_deg$cluster[stat5_deg$logFC>0]),table(stat5_deg$cluster[stat5_deg$logFC<0]))

hits <- findOverlaps(GRanges(stat5_deg[,c(1:4)]),GRanges(tcf1))
stat5_tcf1_deg <- stat5_deg[hits@from,]

stat5_deg$Tcf1_peak <- 0
stat5_deg$Tcf1_peak[hits@from] <- 1

write.xlsx(subset(stat5_deg,logFC>0),'/media/shaoqizhu/easystore/CD8-HP/geneList/Stat5_WT>Stat5_dKO(Tcf1).xlsx')

cbind(table(subset(stat5_tcf1_deg,logFC>0)$cluster),table(subset(stat5_tcf1_deg,logFC<0)$cluster))

cat(as.character(stat5_tcf1_deg$gene_name[stat5_tcf1_deg$cluster==1 & stat5_tcf1_deg$logFC>0]))
