library(ggplot2)
library(ggrepel)
library(stringr)
library(edgeR)
library(grid)
library(pheatmap)

coverage <- read.csv('/media/shaoqizhu/easystore/CD8-HP/H3k27ac/union_peak/coverage.bed',sep='\t',header = F)
merge <- read.csv('/media/shaoqizhu/easystore/CD8-HP/H3k27ac/union_peak/merged.bed',sep='\t',header = F)
name <- read.csv('/media/shaoqizhu/easystore/CD8-HP/H3k27ac/union_peak/name.txt',sep='\t',header = F)
names(coverage)[2:13] <- as.character(name$V1[1:12])
names(merge)<-c('chr','start','end','idx')

RPK <- coverage[,2:13]/coverage[,1]*10^3
TPM <- as.data.frame(t(t(RPK)/apply(RPK, 2, sum)*10^6))
coverage <- coverage[,-1]

pca <- as.data.frame(prcomp(t(TPM[,]),scale.=T)$x)
pca$name <- name$V1[1:12];pca$sample <- name$V1[13:24]

PCA <- prcomp(t(TPM),scale.=T);PC1=(PCA$sdev/sum(PCA$sdev))[1];PC2=(PCA$sdev/sum(PCA$sdev))[2]

ggplot(data = pca,aes(x=PC1, y=PC2, label=name, colour = sample)) + theme_bw() + 
  geom_point() + geom_text_repel(aes(label=name), size = 4, hjust=0, vjust=-1, angle=0) +
  labs(x = paste('PC1 (',round(PC1,3)*100,'% variance)',sep=''),y=paste('PC2 (',round(PC2,3)*100,'% variance)',sep='')) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
library(ggplot2)
library(ggrepel)
library(stringr)
library(edgeR)
library(grid)
library(pheatmap)

coverage <- read.csv('/media/shaoqizhu/easystore/CD8-HP/H3k27ac/union_peak/coverage.bed',sep='\t',header = F)
merge <- read.csv('/media/shaoqizhu/easystore/CD8-HP/H3k27ac/union_peak/merged.bed',sep='\t',header = F)
name <- read.csv('/media/shaoqizhu/easystore/CD8-HP/H3k27ac/union_peak/name.txt',sep='\t',header = F)
names(coverage)[2:13] <- as.character(name$V1[1:12])
names(merge)<-c('chr','start','end','idx')

RPK <- coverage[,2:13]/coverage[,1]*10^3
TPM <- as.data.frame(t(t(RPK)/apply(RPK, 2, sum)*10^6))
coverage <- coverage[,-1]

pca <- as.data.frame(prcomp(t(TPM[,]),scale.=T)$x)
pca$name <- name$V1[1:12];pca$sample <- name$V1[13:24]

PCA <- prcomp(t(TPM),scale.=T);PC1=(PCA$sdev/sum(PCA$sdev))[1];PC2=(PCA$sdev/sum(PCA$sdev))[2]

ggplot(data = pca,aes(x=PC1, y=PC2, label=name, colour = sample)) + theme_bw() + 
  geom_point() + geom_text_repel(aes(label=name), size = 4, hjust=0, vjust=-1, angle=0) +
  labs(x = paste('PC1 (',round(PC1,3)*100,'% variance)',sep=''),y=paste('PC2 (',round(PC2,3)*100,'% variance)',sep='')) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

diff <- numeric(); N = length(unique(pca$exp))
for (i in 1:(N-1)){
  for (j in (i+1):N){
    condition1 = which(pca$exp==unique(pca$exp)[i])
    condition2 = which(pca$exp==unique(pca$exp)[j])
    group <- factor(c(rep(i,length(condition1)),rep(j,length(condition2))))
    y = DGEList(counts = coverage[,c(condition1,condition2)],group = group)
    print(names(coverage[,c(condition1,condition2)]))
    y <- calcNormFactors(y)
    y <- estimateCommonDisp(y)
    fit <- glmQLFit(y)
    qlf <- glmQLFTest(fit)
    fdr <- p.adjust(p = qlf$table$PValue,method = "BH",n = length(qlf$table$PValue))
    diff_idx <- which(fdr < 0.05 & abs(qlf$table$logFC) > log2(1))
    diff <- c(diff,diff_idx)
    change <- cbind(merge[diff_idx,],qlf$table$PValue[diff_idx],fdr[diff_idx],qlf$table$logFC[diff_idx])
    write.table(change,paste('/media/shaoqizhu/easystore/CD8-HP/diff_peaks/',unique(pca$exp)[i],'-',unique(pca$exp)[j],'.bed',sep=''),
                quote = F,sep = '\t',row.names = F,col.names = F)
  }
}
unique(diff)


DNase_coverage <- read.csv('/media/shaoqizhu/easystore/CD8-HP/H3k27ac/DNase_ac/coverage.bed',sep='\t',header = F)
name <- read.csv('/media/shaoqizhu/easystore/CD8-HP/H3k27ac/union_peak/name.txt',sep='\t',header = F)
names(DNase_coverage)[2:13] <- as.character(name$V1[1:12])
swap <- c(2,6,7,8,3,4,5,11,12,13,9,10)
DNase_coverage <- DNase_coverage[,c(1,swap)]
rownames(DNase_coverage) <- rownames(peak_Zscore)

RPK <- DNase_coverage[2:13]/DNase_coverage[,1]*10^3
TPM <- as.data.frame(t(t(RPK)/apply(RPK, 2, sum)*10^6))

H3K27ac_Zscore <- (TPM-apply(X = TPM, MARGIN = 1, mean))/apply(X = TPM, MARGIN = 1, sd)
H3K27ac_Zscore <- H3K27ac_Zscore[rownames(peak_reorder),]

sample_H3K27ac <- data.frame(name$V1[13:24][swap-1]);names(sample_H3K27ac) <- 'condition';
rownames(sample_H3K27ac) <- names(DNase_coverage)[2:13];sample_H3K27ac$condition <- factor(sample_H3K27ac$condition)

pheatmap(H3K27ac_Zscore,cluster_rows=F,cluster_cols = F,show_rownames = F,
         angle_col='315',annotation_row = peak_reorder,annotation_col = sample_H3K27ac,main = 'H3K27ac coverage on DNase diff peak')

DNase_coverage <- read.csv('/media/shaoqizhu/easystore/CD8-HP/H3k27ac/DNase_ac_pooled/coverage.bed',sep='\t',header = F)
rownames(DNase_coverage) <- rownames(peak_Zscore)
name_new <- c('WT_na','dKO_na','WT_s','dKO_s')
names(DNase_coverage)[2:5] <- c('dKO_na','WT_na','dKO_s','WT_s')
DNase_coverage <- DNase_coverage[,c('V1',name_new)]
RPK <- DNase_coverage[2:5]/DNase_coverage[,1]*10^3
TPM <- as.data.frame(t(t(RPK)/apply(RPK, 2, sum)*10^6))

H3K27ac_Zscore <- (TPM-apply(X = TPM, MARGIN = 1, mean))/apply(X = TPM, MARGIN = 1, sd)
H3K27ac_Zscore <- H3K27ac_Zscore[rownames(peak_reorder),]

sample_H3K27ac <- data.frame(name_new);names(sample_H3K27ac) <- 'condition';
rownames(sample_H3K27ac) <- names(DNase_coverage)[2:5];sample_H3K27ac$condition <- factor(sample_H3K27ac$condition)

pheatmap(H3K27ac_Zscore,cluster_rows=F,cluster_cols = F,show_rownames = F,
         angle_col='315',annotation_row = peak_reorder,annotation_col = sample_H3K27ac,main = 'H3K27ac coverage on DNase diff peak')


SE <- read.table('/media/shaoqizhu/easystore/CD8-HP/H3k27ac/bed/coverage.bed')
colnames(SE) <-c("naive_dKO_CD8_201603","naive-dKO-CD8-2_201903","naive-dKO-CD8-2_201908","naive-WT-CD8-1_201903","naive_WT_CD8_201603","naive-WT-CD8-2_201903",
                 "naive-WT-CD8-2_201908","stim-dKO-CD8-2_201903","stim-dKO-CD8-2_201908","stim-WT-CD8-1_201903","stim-WT-CD8-2_201903","stim-WT-CD8-2_201908")
condition <- c("naive-dKO","naive-dKO","naive-dKO","naive-WT","naive-WT","naive-WT",
               "naive-WT","stim-dKO","stim-dKO","stim-WT","stim-WT","stim-WT")

SE_region <- read.table('/media/shaoqizhu/easystore/CD8-HP/H3k27ac/bed/SuperEnhancers_total.bed')

colnames(SE_region) <- names(merge)<-c('chr','start','end')
SE_region$width <- SE_region$end-SE_region$start

plot_pca <- function(coverage,width,condition){ 
  RPK <- coverage/width*10^3
  TPM <- as.data.frame(t(t(RPK)/apply(RPK, 2, sum)*10^6))
  pca <- as.data.frame(prcomp(t(TPM[,]),scale.=T)$x)
  pca$name <- colnames(SE);pca$sample <- condition
  PCA <- prcomp(t(TPM),scale.=T);PC1=(PCA$sdev/sum(PCA$sdev))[1];PC2=(PCA$sdev/sum(PCA$sdev))[2]
  
  ggplot(data = pca,aes(x=PC1, y=PC2, label=name, colour = sample)) + theme_bw() + 
    geom_point() + geom_text_repel(aes(label=name), size = 4, hjust=0, vjust=-1, angle=0) +
    labs(x = paste('PC1 (',round(PC1,3)*100,'% variance)',sep=''),y=paste('PC2 (',round(PC2,3)*100,'% variance)',sep='')) +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
}


diff <- numeric(); N = length(unique(pca$exp))
for (i in 1:(N-1)){
  for (j in (i+1):N){
    condition1 = which(pca$exp==unique(pca$exp)[i])
    condition2 = which(pca$exp==unique(pca$exp)[j])
    group <- factor(c(rep(i,length(condition1)),rep(j,length(condition2))))
    y = DGEList(counts = coverage[,c(condition1,condition2)],group = group)
    print(names(coverage[,c(condition1,condition2)]))
    y <- calcNormFactors(y)
    y <- estimateCommonDisp(y)
    fit <- glmQLFit(y)
    qlf <- glmQLFTest(fit)
    fdr <- p.adjust(p = qlf$table$PValue,method = "BH",n = length(qlf$table$PValue))
    diff_idx <- which(fdr < 0.05 & abs(qlf$table$logFC) > log2(1))
    diff <- c(diff,diff_idx)
    change <- cbind(merge[diff_idx,],qlf$table$PValue[diff_idx],fdr[diff_idx],qlf$table$logFC[diff_idx])
    write.table(change,paste('/media/shaoqizhu/easystore/CD8-HP/diff_peaks/',unique(pca$exp)[i],'-',unique(pca$exp)[j],'.bed',sep=''),
               quote = F,sep = '\t',row.names = F,col.names = F)
  }
}
unique(diff)


DNase_coverage <- read.csv('/media/shaoqizhu/easystore/CD8-HP/H3k27ac/DNase_ac/coverage.bed',sep='\t',header = F)
name <- read.csv('/media/shaoqizhu/easystore/CD8-HP/H3k27ac/union_peak/name.txt',sep='\t',header = F)
names(DNase_coverage)[2:13] <- as.character(name$V1[1:12])
swap <- c(2,6,7,8,3,4,5,11,12,13,9,10)
DNase_coverage <- DNase_coverage[,c(1,swap)]
rownames(DNase_coverage) <- rownames(peak_Zscore)

RPK <- DNase_coverage[2:13]/DNase_coverage[,1]*10^3
TPM <- as.data.frame(t(t(RPK)/apply(RPK, 2, sum)*10^6))

H3K27ac_Zscore <- (TPM-apply(X = TPM, MARGIN = 1, mean))/apply(X = TPM, MARGIN = 1, sd)
H3K27ac_Zscore <- H3K27ac_Zscore[rownames(peak_reorder),]

sample_H3K27ac <- data.frame(name$V1[13:24][swap-1]);names(sample_H3K27ac) <- 'condition';
rownames(sample_H3K27ac) <- names(DNase_coverage)[2:13];sample_H3K27ac$condition <- factor(sample_H3K27ac$condition)

pheatmap(H3K27ac_Zscore,cluster_rows=F,cluster_cols = F,show_rownames = F,
         angle_col='315',annotation_row = peak_reorder,annotation_col = sample_H3K27ac,main = 'H3K27ac coverage on DNase diff peak')

DNase_coverage <- read.csv('/media/shaoqizhu/easystore/CD8-HP/H3k27ac/DNase_ac_pooled/coverage.bed',sep='\t',header = F)
rownames(DNase_coverage) <- rownames(peak_Zscore)
name_new <- c('WT_na','dKO_na','WT_s','dKO_s')
names(DNase_coverage)[2:5] <- c('dKO_na','WT_na','dKO_s','WT_s')
DNase_coverage <- DNase_coverage[,c('V1',name_new)]
RPK <- DNase_coverage[2:5]/DNase_coverage[,1]*10^3
TPM <- as.data.frame(t(t(RPK)/apply(RPK, 2, sum)*10^6))

H3K27ac_Zscore <- (TPM-apply(X = TPM, MARGIN = 1, mean))/apply(X = TPM, MARGIN = 1, sd)
H3K27ac_Zscore <- H3K27ac_Zscore[rownames(peak_reorder),]

sample_H3K27ac <- data.frame(name_new);names(sample_H3K27ac) <- 'condition';
rownames(sample_H3K27ac) <- names(DNase_coverage)[2:5];sample_H3K27ac$condition <- factor(sample_H3K27ac$condition)

pheatmap(H3K27ac_Zscore,cluster_rows=F,cluster_cols = F,show_rownames = F,
         angle_col='315',annotation_row = peak_reorder,annotation_col = sample_H3K27ac,main = 'H3K27ac coverage on DNase diff peak')


SE <- read.table('/media/shaoqizhu/easystore/CD8-HP/H3k27ac/bed/coverage.bed')
colnames(SE) <-c("naive_dKO_CD8_201603","naive-dKO-CD8-2_201903","naive-dKO-CD8-2_201908","naive-WT-CD8-1_201903","naive_WT_CD8_201603","naive-WT-CD8-2_201903",
                 "naive-WT-CD8-2_201908","stim-dKO-CD8-2_201903","stim-dKO-CD8-2_201908","stim-WT-CD8-1_201903","stim-WT-CD8-2_201903","stim-WT-CD8-2_201908")
condition <- c("naive-dKO","naive-dKO","naive-dKO","naive-WT","naive-WT","naive-WT",
               "naive-WT","stim-dKO","stim-dKO","stim-WT","stim-WT","stim-WT")

SE_region <- read.table('/media/shaoqizhu/easystore/CD8-HP/H3k27ac/bed/SuperEnhancers_total.bed')

colnames(SE_region) <- c('chr','start','end')
SE_region$width <- SE_region$end-SE_region$start

plot_pca <- function(coverage,width,condition){ 
  RPK <- coverage/width*10^3
  TPM <- as.data.frame(t(t(RPK)/apply(RPK, 2, sum)*10^6))
  pca <- as.data.frame(prcomp(t(TPM[,]),scale.=T)$x)
  pca$name <- colnames(SE);pca$sample <- condition
  PCA <- prcomp(t(TPM),scale.=T);PC1=(PCA$sdev/sum(PCA$sdev))[1];PC2=(PCA$sdev/sum(PCA$sdev))[2]
  
  ggplot(data = pca,aes(x=PC1, y=PC2, label=name, colour = sample)) + theme_bw() + 
    geom_point() + geom_text_repel(aes(label=name), size = 4, hjust=0, vjust=-1, angle=0) +
    labs(x = paste('PC1 (',round(PC1,3)*100,'% variance)',sep=''),y=paste('PC2 (',round(PC2,3)*100,'% variance)',sep='')) +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
}


plot_pca(SE,SE_region$width,condition)

coverage <- SE
diff <- numeric(); N = length(unique(condition))
for (i in 2){
  for (j in 4){
    condition1 = which(condition==unique(condition)[i])
    condition2 = which(condition==unique(condition)[j])
    group <- factor(c(rep(i,length(condition1)),rep(j,length(condition2))))
    y = DGEList(counts = coverage[,c(condition1,condition2)],group = group)
    print(names(coverage[,c(condition1,condition2)]))
    y <- calcNormFactors(y)
    y <- estimateCommonDisp(y)
    fit <- glmQLFit(y)
    qlf <- glmQLFTest(fit)
    fdr <- p.adjust(p = qlf$table$PValue,method = "BH",n = length(qlf$table$PValue))
    diff_idx <- which(fdr < 0.1 & abs(qlf$table$logFC) > log2(1))
    diff <- c(diff,diff_idx)
    change <- cbind(SE_region[diff_idx,],qlf$table$PValue[diff_idx],fdr[diff_idx],qlf$table$logFC[diff_idx])
    write.table(change,paste('/media/shaoqizhu/easystore/CD8-HP/H3k27ac/diff_SE/',
                             unique(condition)[i],'-',unique(condition)[j],'.bed',sep=''),
                quote = F,sep = '\t',row.names = F,col.names = F)
  }
}
library(GenomicRanges)
hits1 <- findOverlaps(GRanges(change[change[,7]>0,]),GRanges(mm9_deg_50kb))
hits2 <- findOverlaps(GRanges(change[change[,7]<0,]),GRanges(mm9_deg_50kb))
cbind(table(mm9_deg[hits1@to,]$cluster),table(mm9_deg[hits2@to,]$cluster))

cbind(change[change[,7]>0,][hits1@from,1:3][which(mm9_deg[hits1@to,]$cluster==3),],
mm9_deg[hits1@to,][which(mm9_deg[hits1@to,]$cluster==3),])

