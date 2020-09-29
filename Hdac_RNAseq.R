library(ggplot2)
library(ggrepel)
library(stringr)
library(DESeq2)
library(edgeR)
library(grid)
library(factoextra)
library(cluster)
library(GenomicRanges)
library(rtracklayer)
library(pbapply)
library(openxlsx)
library(pheatmap)
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
write.table(exonLen,'/media/shaoqizhu/easystore/Annotation/mm10_exonLen.txt',quote = F,col.names = F,row.names = T)

name <- read.csv('/media/shaoqizhu/easystore/Hdac-RNAseq/info/name.txt',sep='\t',header = F)
coverage <- read.table('/media/shaoqizhu/easystore/Hdac-RNAseq/RNAseq_counts.txt',sep='\t',header = F,row.names=1)
coverage <- coverage[-((nrow(coverage)-4):nrow(coverage)),]
coverage <- coverage[order(rownames(coverage)),]
names(coverage)[1:13] <- as.character(name$V1[1:13])
coverage <- coverage[,c(11:13,3:5,8:9,1:2,10,6:7)]
name$V1[14:26] <- name$V1[14:26][c(11:13,3:5,8:9,1:2,10,6:7)]

coverage <- coverage[,4:13]

RPK <- coverage/exonLen*10^3
RPKM <- as.data.frame(t(t(RPK)/apply(RPK, 2, sum)*10^6))
#write.xlsx(RPKM,paste('/media/shaoqizhu/easystore/Hdac-RNAseq/','Hdac_TPM.xlsx',sep=''),rowNames=T)
RPKM <- RPKM[which(apply(RPKM, 1, var)!=0),]
pca <- as.data.frame(prcomp(t(RPKM),scale.=T)$x)
pca$name <- rownames(pca)
pca$sample <- name$V1[17:26]
PCA <- prcomp(t(RPKM),scale.=T);PC1=(PCA$sdev/sum(PCA$sdev))[1];PC2=(PCA$sdev/sum(PCA$sdev))[2]

ggplot(data = pca,aes(x=PC1, y=PC2, label=name, colour = sample)) + theme_bw() +
  geom_point() + geom_text_repel(aes(label=name), size = 4, hjust=0, vjust=0, angle=0) +
  labs(x = paste('PC1 (',round(PC1,3)*100,'% variance)',sep=''),y=paste('PC2 (',round(PC2,3)*100,'% variance)',sep='')) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())


grep('2',colnames(HdacD01@assays$RNA@counts))

Hdac <- as.data.frame(cbind(D0_WT=apply(HdacD0.data,1,mean),D0_DKO=apply(HdacD0_Dko.data,1,mean),
                            D1_WT=apply(HdacD1.data,1,mean),D1_DKO=apply(HdacD1_Dko.data,1,mean)))
Hdac <- as.data.frame(t(t(Hdac)/apply(Hdac, 2, sum)*10^4))


gene <- intersect(rownames(gene_kmeans),rownames(Hdac))
Hdac_z <- (Hdac[gene,]-apply(X = Hdac[gene,], MARGIN = 1, mean))/apply(X = Hdac[gene,], MARGIN = 1, sd)
Hdac_z <- Hdac_z[Hdac_z$D0_WT!='NaN',]
pheatmap(Hdac_z,show_rownames = T,cluster_rows = F,cluster_cols = F)

plot(log(Hdac[gene,]$D1_WT+0.1),
     log(RPKM[gene,]$`D1-WT1`+0.1))

plot(Hdac$D0_WT,Hdac$D1_DKO)

plot(log(coverage$`D0-WT1`),log(coverage$`D0-DKO1`))

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
    diff_idx <- which(fdr < 0.1 & abs(qlf$table$logFC) > log2(1.5))
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
unique(diff)

DEG <- rownames(coverage)[unique(diff)]



samples <- as.data.frame(pca$sample); names(samples) <- 'sample';rownames(samples) <- pca$name;samples$sample <- factor(samples$sample)

mtx_z <- (RPKM[DEG,]-apply(X = RPKM[DEG,], MARGIN = 1, mean))/apply(X = RPKM[DEG,], MARGIN = 1, sd)
set.seed(2);mtx_kmeans <- kmeans(mtx_z,centers=4,nstart = 25)
gene_kmeans <- as.data.frame(sort(mtx_kmeans$cluster)); names(gene_kmeans) <- 'cluster'; 
gene_kmeans$cluster <- factor(gene_kmeans$cluster)

pheatmap(mtx_z[names(sort(mtx_kmeans$cluster)),],cluster_rows=F,cluster_cols = F,show_rownames = F,
         angle_col='315',annotation_row = gene_kmeans, annotation_col = samples,main = 'Pair-wise DEGs z-score (k-means)')


write.csv(gene_kmeans,'/media/shaoqizhu/easystore/Hdac-RNAseq/gene_kmeans.csv')
  
