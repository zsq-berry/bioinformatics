library(ggplot2)
library(ggrepel)
library(stringr)
library(DESeq2)
library(edgeR)
library(grid)
library(factoextra)
library(cluster)
gtf <- import('/media/shaoqizhu/easystore/LSD1/mm10_2015.gtf')
gtf <- as.data.frame(gtf)
# exon <- gtf[gtf$type=='exon',]
# exonLength <- aggregate(exon$width,list(exon$gene_name),sum)
# exonLength <- exonLength[order(exonLength$Group.1),]
exonLen <- pbsapply(
  X = sort(unique(gtf$gene_name)),FUN = function(x) {
    gene = gtf[gtf$gene_name==x & gtf$type=='exon',]
    gene_ranges <- as.data.frame(reduce(GRanges(gene$seqnames, IRanges(gene$start,gene$end),gene$strand)))
    sum(gene_ranges$width)
  }
)
write.table(exonLen,'/media/shaoqizhu/easystore/LSD1/exonLen.txt',quote = F,col.names = F,row.names = T)

name <- read.csv('/media/shaoqizhu/easystore/LSD1/name.txt',sep='\t',header = F)
coverage <- read.table('/media/shaoqizhu/easystore/LSD1/counts.txt',sep='\t',header = F,row.names=1)
coverage <- coverage[-((nrow(coverage)-4):nrow(coverage)),]
coverage <- coverage[order(rownames(coverage)),]
names(coverage)[1:16] <- as.character(name$V1[1:16])
coverage <- coverage[,-c(8,16)]
coverage <- coverage[,c(5:7,1:4,12:14,8:11)]

RPK <- coverage/exonLen*10^3
RPKM <- as.data.frame(t(t(RPK)/apply(RPK, 2, sum)*10^6))
write.xlsx(RPKM,paste('/media/shaoqizhu/easystore/LSD1/table/','Lsd1_TPM.xlsx',sep=''),rowNames=T)
RPKM <- RPKM[which(apply(RPKM, 1, var)!=0),]
pca <- as.data.frame(prcomp(t(RPKM),scale.=T)$x)
pca$name <- rownames(pca)
pca$sample <- name$V1[c(21:23,17:20,29:31,25:28)]
PCA <- prcomp(t(RPKM),scale.=T);PC1=(PCA$sdev/sum(PCA$sdev))[1];PC2=(PCA$sdev/sum(PCA$sdev))[2]

ggplot(data = pca,aes(x=PC1, y=PC2, label=name, colour = sample)) + theme_bw() +
  geom_point() + #geom_text_repel(aes(label=name), size = 4, hjust=0, vjust=0, angle=0) +
  labs(x = paste('PC1 (',round(PC1,3)*100,'% variance)',sep=''),y=paste('PC2 (',round(PC2,3)*100,'% variance)',sep='')) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

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
    fit <- glmQLFit(y)
    qlf <- glmQLFTest(fit)
    fdr <- p.adjust(p = qlf$table$PValue,method = "BH",n = length(qlf$table$PValue))
    diff_idx <- which(fdr < 0.05 & abs(qlf$table$logFC) > log2(2))
    diff <- c(diff,diff_idx)
    pair_name <- paste(unique(pca$exp)[j],'/',unique(pca$exp)[i],sep='')
    change <- cbind(logFC=qlf$table$logFC,pval=qlf$table$PValue,fdr=fdr)
    colnames(change) <- c(paste(pair_name,'(logFC)'),paste(pair_name,'(pval)'),paste(pair_name,'(fdr)'))
    #write.csv(change,paste('/media/shaoqizhu/easystore/LSD1/table/',unique(pca$exp)[i],'-',unique(pca$exp)[j],'.csv',sep=''),
    #            quote = F,row.names = F)
    #change <- cbind(gene=rownames(coverage)[diff_idx],logFC=qlf$table$logFC[diff_idx],pval=qlf$table$PValue[diff_idx],fdr=fdr[diff_idx])
    # write.table(change,paste('/media/shaoqizhu/easystore/LSD1/DEG/',unique(pca$exp)[i],'-',unique(pca$exp)[j],'.bed',sep=''),
    #             quote = F,sep = '\t',row.names = F,col.names = T)
    # write.table(change[change[,2]>0,],paste('/media/shaoqizhu/easystore/LSD1/DEG/less_than/',unique(pca$exp)[i],'-less-than-',unique(pca$exp)[j],'.bed',sep=''),
    #             quote = F,sep = '\t',row.names = F,col.names = T)
    # write.table(change[change[,2]<0,],paste('/media/shaoqizhu/easystore/LSD1/DEG/more_than/',unique(pca$exp)[i],'-more-than-',unique(pca$exp)[j],'.bed',sep=''),
    #             quote = F,sep = '\t',row.names = F,col.names = T)
  }
}

DEG <- rownames(coverage)[unique(diff)]

samples <- as.data.frame(pca$sample); names(samples) <- 'sample';rownames(samples) <- pca$name;samples$sample <- factor(samples$sample)

mtx_Zscore <- (RPKM[DEG,]-apply(X = RPKM[DEG,], MARGIN = 1, mean))/apply(X = RPKM[DEG,], MARGIN = 1, sd)
mtx_hclust <- hclust(dist(mtx_Zscore), method = "complete")
gene_col <- data.frame(cutree(tree = as.dendrogram(mtx_hclust), k = 4))
names(gene_col)<- 'cluster';gene_col$cluster <- factor(gene_col$cluster)
pheatmap(mtx_Zscore,cluster_rows=T,cluster_cols = F,show_rownames = F,cutree_rows=4,
         angle_col='315',annotation_row = gene_col, annotation_col = samples,main = 'Pair-wise DEGs z-score (hclust)')

set.seed(2);

fviz_nbclust(mtx_Zscore, kmeans, method = "wss")

set.seed(2);mtx_kmeans <- kmeans(mtx_Zscore,centers=4,nstart = 25)
gene_kmeans <- as.data.frame(sort(mtx_kmeans$cluster)); names(gene_kmeans) <- 'cluster'; 
gene_kmeans$cluster <- factor(gene_kmeans$cluster)

pheatmap(mtx_Zscore[names(sort(mtx_kmeans$cluster)),],cluster_rows=F,cluster_cols = F,show_rownames = F,
         angle_col='315',annotation_row = gene_kmeans, annotation_col = samples,main = 'Pair-wise DEGs z-score (k-means)')

write.csv(gene_kmeans,'/media/shaoqizhu/easystore/LSD1/DEG/gene_kmeans.csv')

gene_kmeans$gene_name <- rownames(gene_kmeans)
gene_kmeans <- gene_kmeans[order(gene_kmeans$gene_name),]
TPM <- read.xlsx('/media/shaoqizhu/easystore/LSD1/table/Lsd1_TPM.xlsx',rowNames = T)
TPM <- TPM[rownames(TPM)%in%DEG,]
TPM$kmeans_cluster <- gene_kmeans$cluster
TPM <- TPM[,c(1:14,33,15:32)]
TPM <- TPM[order(TPM$kmeans_cluster),]
write.xlsx(TPM,'/media/shaoqizhu/easystore/LSD1/table/Lsd1_DEG_TPM.xlsx',rowNames = T)

overlap <- matrix(data=0,nrow = 5,ncol = 5)
for (i in 1:5){
  for (j in 1:5){
    overlap[i,j] <- length(intersect(rownames(gene_kmeans)[gene_kmeans$cluster==i],
                                     rownames(gene_col)[gene_col$cluster==j]))
  }
}
overlap <- as.data.frame(overlap)
names(overlap) <- 1:5
pheat <- pheatmap(overlap,display_numbers = TRUE,cluster_cols = F,cluster_rows = F,
                  angle_col=0,treeheight_row=10,treeheight_col=10,number_format = "%.0f",fontsize_number = 12,fontsize = 12)
table(gene_kmeans$cluster)
table(gene_col$cluster)

{File = '4'
  GO <- read.csv(paste('/media/shaoqizhu/easystore/LSD1/DAVID/',File,'.txt',sep=''),header = T, sep = "\t")
  GO$Term <- str_split_fixed(GO$Term,'~',2)[,2]
  GO$Term <- factor(GO$Term, levels=rev(unique(GO$Term)))
  GO <-GO[1:20,]
  ggplot(data = GO, aes(x = GO$Count,y = GO$Term),) +
    geom_point(size = 3,aes(color = -log10(GO$FDR))) +
    scale_color_gradient(low="green",high="red") +
    labs(color = expression(-log[10](FDR)),x = "Gene Count",y="GO term",
         title=paste('DAVID GO in kmeans cluster',File)) +
    theme(plot.title=element_text(size=12,vjust=0.5,hjust=0.5),
          legend.title=element_text(size=10),
          legend.text=element_text(size=10),
          legend.position="right",
          axis.text=element_text(size=10),
          axis.title=element_text(size=10),)
}
