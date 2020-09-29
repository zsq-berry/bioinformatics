library(ggplot2)
library(ggrepel)
library(stringr)
library(DESeq2)
library(edgeR)
library(grid)

coverage <- read.csv('/media/shaoqizhu/easystore/CD8-HP/counts/coverage.bed',sep='\t',header = F)
merge <- read.csv('/media/shaoqizhu/easystore/CD8-HP/counts/merged.bed',sep='\t',header = F)
name <- read.csv('/media/shaoqizhu/easystore/CD8-HP/counts/name.txt',sep='\t',header = F)
name$V1 <- gsub('-','_',name$V1)
names(coverage)[2:11] <- as.character(name$V1[1:10])
names(merge)<-c('chr','start','end','idx')
coverage$V1 <- merge$end-merge$start

RPK <- coverage[,2:11]/coverage[,1]*10^3
RPKM <- as.data.frame(t(t(RPK)/apply(RPK, 2, sum)*10^6))
coverage <- coverage[,-1]

pca <- as.data.frame(prcomp(t(RPKM[,]),scale.=T)$x)
pca$name <- name$V1[1:10];pca$exp <- name$V1[11:20]

ggplot(data = pca,aes(x=PC1, y=PC2, label=name, colour = exp)) + 
  geom_point() + geom_text_repel(aes(label=name), size = 4, hjust=0, vjust=-1, angle=0)


PCA <- prcomp(t(coverage[,]),scale.=T)
PCA <- prcomp(t(RPKM[,]),scale.=T)
(PCA$sdev/sum(PCA$sdev))[1:2]

diff <- numeric(); N = length(unique(pca$exp));#Name <- character()
for (i in 1:(N-1)){
  for (j in (i+1):N){
    condition1 = which(pca$exp==unique(pca$exp)[i])
    condition2 = which(pca$exp==unique(pca$exp)[j])
    coldata <- as.data.frame(pca$exp[c(condition1,condition2)])
    colnames(coldata) = 'condition'
    rownames(coldata) = as.character(names(coverage)[c(condition1,condition2)])
    dds <- DESeqDataSetFromMatrix(countData = coverage[,c(condition1,condition2)],
                                  colData = coldata,
                                  design = ~ condition)
    dds <- DESeq(dds,fitType = 'local')
    res <- results(dds)
    print(names(coverage[,c(condition1,condition2)]))
    #diff=rbind(diff,(length(which(res$pvalue < 0.01 & abs(res$log2FoldChange) > log2(3)))))
    diff_idx <- which(res$padj < 0.01 & abs(res$log2FoldChange) > log2(3))
    diff <- c(diff,diff_idx)
    #change = cbind(res$log2FoldChange,res$padj)
    #Name <- c(Name,paste(unique(pca$exp)[j],'/',unique(pca$exp)[i],'(log2FC)',sep=''),paste(unique(pca$exp)[j],'/',unique(pca$exp)[i],'(padj)',sep=''))
    #diff <- c(diff,diff_0)
    change <- cbind(merge[diff_idx,],res$log2FoldChange[diff_idx])
    #write.table(change,paste('/media/shaoqizhu/easystore/hdac/merge_cut_150/diff_peaks/',unique(pca$exp)[i],'-',unique(pca$exp)[j],'.bed',sep=''),
    #            quote = F,sep = '\t',row.names = F,col.names = F)
  }
}

length(unique(diff))

diff <- as.data.frame(diff)
names(diff) <- Name; names(merge) <- c('chrome','start','end','idx');names(RPKM) <- name$V1
write.table(cbind(merge,round(RPKM,2),diff[,c(3,4,5,6,15,16)]),paste('/media/shaoqizhu/easystore/hdac/bampe/','merge_peak_diff.bed',sep=''),
            quote = F,sep = '\t',row.names = F,col.names = T)

write.table(merge[unique(diff),],paste('/media/shaoqizhu/easystore/hdac/pca_no_dkmyc/','merge_diff_peak','.bed',sep=''),
            quote = F,sep = '\t',row.names = F,col.names = F)

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
    diff_idx <- which(fdr < 0.05 & abs(qlf$table$logFC) > log2(2))
    diff <- c(diff,diff_idx)
    change <- cbind(merge[diff_idx,],qlf$table$PValue[diff_idx],fdr[diff_idx],qlf$table$logFC[diff_idx])
    #write.table(change,paste('/media/shaoqizhu/easystore/CD8-HP/diff_peaks/',unique(pca$exp)[i],'-',unique(pca$exp)[j],'.bed',sep=''),
    #           quote = F,sep = '\t',row.names = F,col.names = F)
  }
}
unique(diff)


union_diff <- read.csv('/media/shaoqizhu/easystore/CD8-HP/diff_peaks/union_diff_peaks.txt',sep='\t',header = F)
union_diff_rpkm <- RPKM[union_diff$V1,]
peak_Zscore <- (union_diff_rpkm-apply(X = union_diff_rpkm, MARGIN = 1, mean))/apply(X = union_diff_rpkm, MARGIN = 1, sd)

rename <- c("WT_na1","WT_na2","11_WT_na","12_dKO_na","dKO_na1","5_WT_s","WT_s2","10_dKO_s","8_dKO_s","dKO_s1",
            "WT_na","WT_na","WT_na","dKO_na","dKO_na","WT_s","WT_s","dKO_s","dKO_s","dKO_s")

peak_Zscore <- peak_Zscore[,rename[1:10]]

peak_hclust <- hclust(dist(peak_Zscore), method = "complete")
peak_col <- data.frame(cutree(tree = as.dendrogram(peak_hclust), k = 7))

names(peak_col)<- 'cluster';peak_col$cluster <- factor(peak_col$cluster)

sample_DNase <- data.frame(rename[11:20]);names(sample_DNase) <- 'condition';
rownames(sample_DNase) <- rename[1:10];sample_DNase$condition <- factor(sample_DNase$condition)

pheatmap(peak_Zscore,cluster_rows=T,cluster_cols = F,show_rownames = F,cutree_rows=7,
         angle_col='315',annotation_row = peak_col, annotation_col = sample_DNase,main = 'Pair-wise diff peaks z-score (hclust)')

#fviz_nbclust(peak_Zscore, kmeans, method = "wss")

set.seed(2);matrix_kmeans <- kmeans(peak_Zscore,centers=7,nstart=25)
peak_kmeans <- as.data.frame(sort(matrix_kmeans$cluster)); names(peak_kmeans) <- 'cluster'; 

reorder_peak <- c(4,3,6,7,2,5,1);peak_reorder <- peak_kmeans;
for (i in 1:7) {peak_reorder$cluster[peak_kmeans$cluster==reorder_peak[i]] <- i}
name_ordered <- rownames(peak_kmeans)[order(peak_reorder$cluster)]
peak_reorder <- as.data.frame(sort(peak_reorder$cluster));names(peak_reorder)<- 'cluster'
rownames(peak_reorder) <- name_ordered
peak_reorder$cluster <- factor(peak_reorder$cluster)

pheatmap(peak_Zscore[rownames(peak_reorder),],cluster_rows=F,cluster_cols = F,show_rownames = F,
         angle_col='315',annotation_row = peak_reorder, annotation_col = sample_DNase,main = 'Pair-wise diff peak z-score (k-means)')

peak_diff_ac <- read.table('/media/shaoqizhu/easystore/CD8-HP/H3k27ac/DNase_diff/coverage.bed')
rownames(peak_diff_ac) <- rownames(peak_Zscore)
pheatmap(log(peak_diff_ac[rownames(peak_reorder),c(2:12)]+1),cluster_rows=F,cluster_cols = F,show_rownames = F)
         
overlap <- matrix(data=0,nrow = 7,ncol = 7)
for (i in 1:7){
  for (j in 1:7){
    overlap[i,j] <- length(intersect(rownames(peak_kmeans)[peak_kmeans$cluster==i],
                                     rownames(peak_col)[peak_col$cluster==j]))
  }
}
overlap <- as.data.frame(overlap)
names(overlap) <- 1:7
pheat <- pheatmap(overlap,display_numbers = TRUE,cluster_rows = F,cluster_cols = F,
                  angle_col=0,treeheight_row=10,treeheight_col=10,number_format = "%.0f",fontsize_number = 12,fontsize = 12)

table(peak_kmeans$cluster)
table(peak_col$cluster)

peak_idx <- matrix_kmeans$cluster
for (i in 1:7) {peak_idx[matrix_kmeans$cluster==reorder_peak[i]] <- i}

diff_peak_kmeans <- cbind(merge[union_diff$V1,],peak_idx)
write.table(cbind(merge[union_diff$V1,],cluster_idx),'/media/shaoqizhu/easystore/CD8-HP/bedFile/diff_peak_kmeans.bed',
            quote = F,row.names = F,col.names = F,sep = '\t')
#bedtools intersect -wo -a diff_peak_kmeans.bed -b mm9_2014_reduced_50kb.bed | cut -f1,2,3,4,5,10 > diff_peak_kmeans_genes.bed
write.table(diff_peak_kmeans[diff_peak_kmeans$cluster_idx==1,1:4],'/media/shaoqizhu/easystore/CD8-HP/geneList/diff_peak_cluster1.bed',
            quote = F,row.names = F,col.names = F,sep = '\t')
as.vector(GO[2,])
{
  File = 'diff_peak_cluster1'
  GO <- read.csv(paste('/media/shaoqizhu/easystore/CD8-HP/GREAT/',File,'.tsv',sep=''),header = T, sep = "\t")
  colnames(GO)[1] <- "Term"
  GO$Term <- factor(GO$Term, levels=rev(unique(GO$Term)))
  ggplot(data = GO, aes(x = -log10(Binom.FDR.Q.Val),y = Term),) +
    geom_point(aes(colour = Binom.Fold.Enrichment,size = Binom.Observed.Region.Hits)) +
    scale_color_gradient(low="green",high="red") +
    labs(colour = 'Fold.Enrichment',size = 'Observed.Region.Hits',
         x = "FDR",y="GO term",
         title=paste('GREAT GO in cluster','1')) +
    theme(plot.title=element_text(size=12,vjust=0.5,hjust=0.5),
          legend.title=element_text(size=10),
          legend.text=element_text(size=10),
          legend.position="right",
          axis.text=element_text(size=10),
          axis.title=element_text(size=10),)
}

hits <- suppressWarnings(findOverlaps(GRanges(diff_peak_kmeans),GRanges(mm9_deg_50kb)))
gene_overlap <- cbind(diff_peak_kmeans[hits@from,],mm9_deg_50kb[hits@to,])
gene_overlap <- as.data.frame(cbind(gene_overlap$cluster,gene_overlap$cluster_idx))
pheatmap(table(gene_overlap),display_numbers = TRUE,number_format = "%.0f",fontsize_number = 12,
         cluster_rows = F,cluster_cols =F,angle_col=0,fontsize = 13)
write.csv(gene_overlap,'/media/shaoqizhu/easystore/CD8-HP/geneList/gene_overlap_DNase_RNAseq.csv',row.names = F)

peak_kmeans_gene <- read.table('/media/shaoqizhu/easystore/CD8-HP/geneList/diff_peak_kmeans_genes.bed')
#unique_gene <- names(table(peak_kmeans_gene$V6))[which(table(peak_kmeans_gene$V6)==1)]
#peak_kmeans_gene <- peak_kmeans_gene[peak_kmeans_gene$V6 %in% unique_gene,]
overlap_gene <- matrix(data=0,nrow=8,ncol=7); gene_name <- as.character(); gene_idx <- matrix(NA,nrow=0,ncol=2)
for (i in 1:8){
  for (j in 1:7){
    overlap_gene[i,j] <- length(intersect(rownames(gene_reorder)[gene_reorder$cluster==i],
                              peak_kmeans_gene$V6[peak_kmeans_gene$V5==j]))
    temp <- intersect(rownames(gene_reorder)[gene_reorder$cluster==i],
              peak_kmeans_gene$V6[peak_kmeans_gene$V5==j])
    gene_idx <- rbind(gene_idx,cbind(rep(i,length(temp)),rep(j,length(temp))))
    gene_name <- c(gene_name,temp)
  }
}

overlap_gene <- as.data.frame(overlap_gene)
names(overlap_gene) <- 1:7
pheatmap(overlap_gene,display_numbers = TRUE,number_format = "%.0f",fontsize_number = 12,
         cluster_rows = F,cluster_cols =F,angle_col=0,fontsize = 13)
table(gene_reorder$cluster)
table(peak_kmeans_gene$V5)
overlap_gene/norm

norm <- apply(overlap_gene,1,sum)%*%t(apply(overlap_gene,2,sum))/sum(overlap_gene)

gene_list <- cbind(gene_name,gene_idx);colnames(gene_list)[2:3] <- c('RNAseq_cluster','DNase_cluster')
write.csv(gene_list,'/media/shaoqizhu/easystore/CD8-HP/geneList/gene_overlap_DNase_RNAseq.csv',row.names = F)
gene_name[1]


peak_mean <- t(sapply(X=1:7, function(x){apply(peak_Zscore[matrix_kmeans$cluster==reorder_peak[x],],2,mean)}))
deg_mean <- t(sapply(X=1:8, function(x){apply(mtx_Zscore[mtx_kmeans$cluster==reorder_gene[x],],2,mean)}))
deg_mean <- deg_mean[,-c(6,7)]
correlation <- as.data.frame(cor(t(deg_mean),t(peak_mean)))
colnames(correlation) <- 1:7; rownames(correlation) <- 1:8
pheatmap(correlation,display_numbers = TRUE,number_format = "%.2f",fontsize_number = 12,
         cluster_rows = F,cluster_cols =F,angle_col=0,fontsize = 13)

cor_para <- round(cor(as.vector(t(correlation)),c(t(overlap_gene/norm))),4)
ggplot()+geom_point(aes(x=as.vector(t(correlation)),y=c(t(overlap_gene/norm))))+
  labs(x='correlation',y='gene overlap',title = paste('cor =',cor_para))+
  theme(plot.title=element_text(size=12,vjust=0.5,hjust=0.5))

cor_para <- round(cor(c(t(table(peak_overlap))),c(t(overlap_gene))),4)
ggplot()+geom_point(aes(x=c(t(table(peak_overlap))),y=c(t(overlap_gene))))+
  labs(x='peak overlap',y='gene overlap',title = paste('cor =',cor_para))+
  theme(plot.title=element_text(size=12,vjust=0.5,hjust=0.5))

set.seed(1.5);matrix_kmeans <- kmeans(peak_Zscore,centers=7)
set.seed(5);matrix_kmeans2 <- kmeans(peak_Zscore,centers=7)
overlap0 <- matrix(data=0,nrow = 7,ncol = 7)
for (i in 1:7){
  for (j in 1:7){
    overlap0[i,j] <- length(intersect(names(matrix_kmeans$cluster)[matrix_kmeans$cluster==i],
                                      names(matrix_kmeans2$cluster)[matrix_kmeans2$cluster==j]))
  }
}
overlap0 <- as.data.frame(overlap0);names(overlap0) <- 1:7
pheatmap(overlap0,display_numbers = TRUE,cluster_rows = T,cluster_cols = T,
         angle_col=0,treeheight_row=10,treeheight_col=10,number_format = "%.0f",fontsize_number = 12,fontsize = 12)


diff_peak_kmeans_summit <- read.table('/media/shaoqizhu/easystore/CD8-HP/macs2/summits/diff_peak_kmeans_summit.bed')

summit_kmeans <- character()
for (i in unique(diff_peak_kmeans_summit$V4)){
  summit <- diff_peak_kmeans_summit[diff_peak_kmeans_summit$V4 == i,c(7,8,9,10,5,11)]
  summit_max <- summit[summit$V11==max(summit$V11),c(1,2,3,4,5)]
  summit_kmeans <- rbind(summit_kmeans,summit_max)
}

for (i in 1:7){
  write.table(summit_kmeans[summit_kmeans$V5==i,1:4],
              paste('/media/shaoqizhu/easystore/CD8-HP/macs2/summit_cluster',i,'.bed',sep = ''),
              quote = F,col.names = F,row.names = F,sep = '\t')
}
summit_max

hits <- findOverlaps(GRanges(ctcf),GRanges(diff_peak_kmeans))

hits1 <- findOverlaps(GRanges(tcf1),GRanges(diff_peak_kmeans[hits@to,]))
cbind(table(diff_peak_kmeans$cluster_idx[hits@to]),
      table(diff_peak_kmeans$cluster_idx[hits@to][hits1@to]))

tcf1 <- read.table('/media/shaoqizhu/easystore/CD8-HP/CD8_HP_naive/NaiveCD8_TCF1_BAMPE_peaks_filtered.bed')
stat5_dko_24h <- read.table('/media/shaoqizhu/easystore/Stat5/narrowPeak/dko_24h_1_peaks.narrowPeak')
stat5_wt_24h <- read.table('/media/shaoqizhu/easystore/Stat5/narrowPeak/wt_24h_1_peaks.narrowPeak')
names(tcf1)[1:3] <- c('chr','start','end');
names(stat5_dko_24h)[1:3] <- c('chr','start','end');
names(stat5_wt_24h)[1:3] <- c('chr','start','end');
tcf1_DNase <- findOverlaps(GRanges(tcf1),GRanges(diff_peak_kmeans))
stat5_dko_DNase <- findOverlaps(GRanges(stat5_dko_24h),GRanges(diff_peak_kmeans))
stat5_wt_DNase <- findOverlaps(GRanges(stat5_wt_24h),GRanges(diff_peak_kmeans))

cbind(table(diff_peak_kmeans$cluster_idx[tcf1_DNase@to]),
      table(diff_peak_kmeans$cluster_idx[stat5_wt_DNase@to]),
      table(diff_peak_kmeans$cluster_idx[stat5_dko_DNase@to]))

cbind(table(diff_peak_kmeans$cluster_idx[intersect(tcf1_DNase@to,stat5_wt_DNase@to)]),
      table(diff_peak_kmeans$cluster_idx[intersect(tcf1_DNase@to,stat5_dko_DNase@to)]))

setdiff(stat5_dko_DNase@from,stat5_wt_DNase@from)

overlap_stat5_dko_DNase <- cbind(stat5_dko_24h[stat5_dko_DNase@from,c(1:4)],diff_peak_kmeans[stat5_dko_DNase@to,c(4,5)])
overlap_stat5_wt_DNase <- cbind(stat5_wt_24h[stat5_wt_DNase@from,c(1:4)],diff_peak_kmeans[stat5_wt_DNase@to,c(4,5)])

hits1 <- findOverlaps(GRanges(overlap_stat5_dko_DNase),GRanges(tcf1))
hits2 <- findOverlaps(GRanges(overlap_stat5_wt_DNase),GRanges(tcf1))

overlap_stat5_dko_DNase_tcf1 <- overlap_stat5_dko_DNase[hits1@from,]
overlap_stat5_wt_DNase_tcf1 <- overlap_stat5_wt_DNase[hits2@from,]
cbind(table(overlap_stat5_wt_DNase_tcf1$cluster_idx),
      table(overlap_stat5_dko_DNase_tcf1$cluster_idx))

write.table(overlap_stat5_wt_DNase_tcf1,'/media/shaoqizhu/easystore/CD8-HP/bedFile/overlap_stat5_wt_DNase_tcf1.bed',
            quote = F,col.names = F,row.names = F,sep = '\t')

for (i in 7){
  gr1 <- overlap_stat5_wt_DNase_tcf1[overlap_stat5_wt_DNase_tcf1$cluster_idx==i,]
  gr2 <- overlap_stat5_dko_DNase_tcf1[overlap_stat5_dko_DNase_tcf1$cluster_idx==i,]
  hits <- findOverlaps(GRanges(gr1),GRanges(gr2))
  print(c(nrow(gr1[-hits@from,]),length(hits@from),nrow(gr2[-hits@to,])))
  hits2 <- suppressWarnings(findOverlaps(GRanges(gr1[-hits@from,]),GRanges(mm9_deg_50kb)))
  gr1_deg <- cbind(gr1[-hits@from,][hits2@from,],mm9_deg[hits2@to,])
  table(gr1_deg$cluster)
}
cat(as.character(gr1_deg$gene_name[gr1_deg$cluster==1]))


hits <- findOverlaps(GRanges(diff_peak_kmeans),GRanges(mm9_promoter))
diff_peak_cluster1 <- diff_peak_kmeans[hits@from,][diff_peak_kmeans[hits@from,]$cluster_idx==1,]
diff_peak_cluster5 <- diff_peak_kmeans[hits@from,][diff_peak_kmeans[hits@from,]$cluster_idx==5,]

mm9_size <- read.table('/media/shaoqizhu/easystore/read-through/annotation/mm9.genome')
mm9_size <- mm9_size[1:22,]
rownames(mm9_size) <- mm9_size$V1
mm9_size <- mm9_size[order(mm9_size$V1),]

write.table(cbind(gsub(' ','',paste(str_split_fixed(diff_peak_cluster1$chr,'r',2)[,2],
                  format(floor(diff_peak_cluster1$start/10000)*10000, scientific=F),
                  format(floor(diff_peak_cluster1$start/10000)*10000, scientific=F),sep=':')),
                  paste(str_split_fixed(mm9_size[diff_peak_cluster1$chr,]$V1,'r',2)[,2],0,
                        mm9_size[diff_peak_cluster1$chr,]$V2,sep = ':')),'/media/shaoqizhu/easystore/HiC_CD8HP/cluster1_prom_region.txt'
            ,sep = ' ',quote = F,col.names = F,row.names = F)

#i=0; while IFS= read -r line; do let "i=i+1"; java -jar juicer_tools.jar dump observed KR WT_CD8_1910_Juicebox.hic $line BP 10000 ${i}.txt; done < cluster1_prom_region.txt
#i=0; while IFS= read -r line; do let "i=i+1"; java -jar juicer_tools.jar dump observed KR DKO_CD8_Juicebox_2016.hic $line BP 10000 DKO_2016/$(cut -d':' -f 1 <<<$line).txt; done < mm9_size.txt

write.table(cbind(paste(str_split_fixed(mm9_size$V1,'r',2)[,2],0,mm9_size$V2,sep = ':'),
                  paste(str_split_fixed(mm9_size$V1,'r',2)[,2],0,mm9_size$V2,sep = ':')),
            '/media/shaoqizhu/easystore/HiC_CD8HP/mm9_size.txt',
            sep = ' ',quote = F,col.names = F,row.names = F)

diff_peak_1 <- cbind(str_split_fixed(diff_peak_cluster1$chr,'r',2)[,2],
      gsub(' ','',format(floor(diff_peak_cluster1$start/10000)*10000, scientific=F)))
write.table(diff_peak_1,'/media/shaoqizhu/easystore/HiC_CD8HP/cluster1_prom_region.txt',sep = ' ',quote = F,col.names = F,row.names = F)

diff_peak_5 <- cbind(str_split_fixed(diff_peak_cluster5$chr,'r',2)[,2],
                     gsub(' ','',format(floor(diff_peak_cluster5$start/10000)*10000, scientific=F)))
write.table(diff_peak_5,'/media/shaoqizhu/easystore/HiC_CD8HP/cluster5_prom_region.txt',sep = ' ',quote = F,col.names = F,row.names = F)

deg_1 <- cbind(str_split_fixed(mm9_deg_prom$chr[mm9_deg_prom$cluster==1],'r',2)[,2],
               gsub(' ','',format(floor((mm9_deg_prom$start[mm9_deg_prom$cluster==1]+1000)/10000)*10000, scientific=F)))
write.table(deg_1,'/media/shaoqizhu/easystore/HiC_CD8HP/cluster1_deg_prom.txt',sep = ' ',quote = F,col.names = F,row.names = F)

deg_5 <- cbind(str_split_fixed(mm9_deg_prom$chr[mm9_deg_prom$cluster==5],'r',2)[,2],
               gsub(' ','',format(floor((mm9_deg_prom$start[mm9_deg_prom$cluster==5]+1000)/10000)*10000, scientific=F)))
write.table(deg_5,'/media/shaoqizhu/easystore/HiC_CD8HP/cluster5_deg_prom.txt',sep = ' ',quote = F,col.names = F,row.names = F)

deg_6 <- cbind(str_split_fixed(mm9_deg_prom$chr[mm9_deg_prom$cluster==6],'r',2)[,2],
               gsub(' ','',format(floor((mm9_deg_prom$start[mm9_deg_prom$cluster==6]+1000)/10000)*10000, scientific=F)))
write.table(deg_6,'/media/shaoqizhu/easystore/HiC_CD8HP/cluster6_deg_prom.txt',sep = ' ',quote = F,col.names = F,row.names = F)

mm9_promoter[mm9_promoter$gene_name%in%WT_na_lt_DKO_na,]
mm9_promoter[mm9_promoter$gene_name%in%WT_na_mt_DKO_na,]
write.table(cbind(str_split_fixed(mm9_promoter[mm9_promoter$gene_name%in%WT_na_mt_DKO_na,1],'r',2)[,2],
      gsub(' ','',format(floor((mm9_promoter[mm9_promoter$gene_name%in%WT_na_mt_DKO_na,2]+1000)/10000)*10000, scientific=F))),
      '/media/shaoqizhu/easystore/HiC_CD8HP/na_WT_mt_DKO_prom.txt',sep = ' ',quote = F,col.names = F,row.names = F)
write.table(cbind(str_split_fixed(mm9_promoter[mm9_promoter$gene_name%in%WT_na_lt_DKO_na,1],'r',2)[,2],
      gsub(' ','',format(floor((mm9_promoter[mm9_promoter$gene_name%in%WT_na_lt_DKO_na,2]+1000)/10000)*10000, scientific=F))),
      '/media/shaoqizhu/easystore/HiC_CD8HP/na_WT_lt_DKO_prom.txt',sep = ' ',quote = F,col.names = F,row.names = F)

#cat WT_2016/1.txt | grep $(head -1 cluster1_prom_region.txt | cut -d' ' -f2) | awk '{if($2-$1>10000 || $1-$2>10000) print$0}' | awk -F '\t' '{sum += $3} END {print sum}'
#echo 'sum' > WT_2016_sum.txt;while IFS= read -r line; do i=$(cat WT_2016/"$(cut -d' ' -f1 <<<$line)".txt | grep $(cut -d' ' -f2 <<<$line) | awk '{if($2-$1>10000 || $1-$2>10000) print$0}' | grep -v NaN | awk -F '\t' '{sum += $3} END {print sum}'); sed -i "$ a "$i"" WT_2016_sum.txt;; done < cluster1_prom_region.txt

hit_sum <- cbind(WT=read.table('/media/shaoqizhu/easystore/HiC_CD8HP/out/naive_DKO_2016_sum_cluster1_deg_prom.txt',header = T)[,1],
                  DKO=read.table('/media/shaoqizhu/easystore/HiC_CD8HP/out/stim_DKO_2016_sum_cluster1_deg_prom.txt',header = T)[,1])
ratio <- (hit_sum[,2]-hit_sum[,1])/hit_sum[,1]
hist(ratio,breaks=seq(min(ratio),max(ratio),length.out=20),
     main='Sum of contact density DKO (stim-naive)/naive for deg cluster1',xlab='(stim-naive)/naive')
boxplot(ratio,main='Sum of contact density (DKO-WT)/WT cluster 5 deg',xlab='(DKO-WT)/WT')
    