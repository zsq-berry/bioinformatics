coverage <- read.csv('/media/shaoqizhu/easystore/CD8-HP/CTCF/union_peak/no_24h/coverage.bed',sep='\t',header = F)
merge <- read.csv('/media/shaoqizhu/easystore/CD8-HP/CTCF/union_peak/no_24h/merged.bed',sep='\t',header = F)
name <- c("DKO_72h_1","DKO_72h_2","DKO_na_1","DKO_na_2",
          "WT_72h_1","WT_72h_2","WT_na_1","WT_na_2")
name1 <- c("DKO_72h","DKO_72h","DKO_na","DKO_na",
          "WT_72h","WT_72h","WT_na","WT_na")

names(coverage)[-1] <- as.character(name)
names(merge)<-c('chr','start','end','idx')

RPK <- coverage[,-1]/coverage[,1]*10^3
RPKM <- as.data.frame(t(t(RPK)/apply(RPK, 2, sum)*10^6))
coverage <- coverage[,-1]

pca <- as.data.frame(prcomp(t(RPKM[,]),scale.=T)$x)
pca$name <- colnames(coverage);pca$exp <- name1

ggplot(data = pca,aes(x=PC1, y=PC2, label=name, colour = exp)) + 
  geom_point() + geom_text_repel(aes(label=name), size = 4, hjust=0, vjust=-1, angle=0)

PCA <- prcomp(t(RPKM[,]),scale.=T)
(PCA$sdev/sum(PCA$sdev))[1:2]

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
    write.table(change,paste('/media/shaoqizhu/easystore/CD8-HP/CTCF/peak_diff_all/pair-wise-no_24h/',unique(pca$exp)[i],'-',unique(pca$exp)[j],'.bed',sep=''),
               quote = F,sep = '\t',row.names = F,col.names = F)
  }
}
length(unique(diff))
write.table(merge[unique(diff),],
            paste('/media/shaoqizhu/easystore/CD8-HP/CTCF/peak_diff_all/merged_diff_peaks.bed',sep = ''),
            quote = F,row.names = F,col.names = F,sep = '\t')

peak_diff <- RPKM[unique(diff),c(7,8,3,4,5,6,1,2)]
peak_diff_z <- (peak_diff-apply(X = peak_diff, MARGIN = 1, mean))/apply(X = peak_diff, MARGIN = 1, sd)

set.seed(2);peak_diff_kmeans <- kmeans(peak_diff_z,centers=6,nstart = 25)
peak_kmeans <- as.data.frame(sort(peak_diff_kmeans$cluster)); names(peak_kmeans) <- 'cluster'; 
peak_kmeans$cluster <- factor(peak_kmeans$cluster)
pheatmap::pheatmap(peak_diff_z[rownames(peak_kmeans),],cluster_rows=F,cluster_cols = F,show_rownames = F,
                   angle_col='315',annotation_row = peak_kmeans, main = 'RNA-seq z-score (k-means)')

reorder_peak <- c(3,1,4,2,6,5);peak_reorder <- peak_kmeans;
for (i in 1:6) {peak_reorder$cluster[peak_kmeans$cluster==reorder_peak[i]] <- i}
name_ordered <- rownames(peak_kmeans)[order(peak_reorder$cluster)]
peak_reorder <- as.data.frame(sort(peak_reorder$cluster));names(peak_reorder)<- 'cluster'
rownames(peak_reorder) <- name_ordered
peak_reorder$cluster <- factor(peak_reorder$cluster)
pheatmap::pheatmap(peak_diff_z[rownames(peak_reorder),],cluster_rows=F,cluster_cols = F,show_rownames = F,
                   angle_col='315',annotation_row = peak_reorder,main = 'CTCF diff peak z-score (k-means)')

for (i in 1:6){
  hits <- findOverlaps(GRanges(Tcf1),GRanges(merge[rownames(peak_reorder)[peak_reorder$cluster==i],]))
  #print(c(length(unique(hits@from)),hits@nRnode))
  print(hits@nRnode)
}

for (i in 1:6){
  write.table(merge[rownames(peak_reorder)[peak_reorder$cluster==i],],
              paste('/media/shaoqizhu/easystore/CD8-HP/CTCF/peak_diff_all//diff_kmeans_no24h/cluster_',i,'.bed',sep = ''),
              quote = F,row.names = F,col.names = F,sep = '\t')
}

CTCF_cluster <- merge[rownames(peak_reorder)[peak_reorder$cluster==1],]
hits <- findOverlaps(GRanges(Tcf1),GRanges(CTCF_cluster))
write.table(CTCF_cluster[-unique(hits@to),],
            '/media/shaoqizhu/easystore/CD8-HP/CTCF/peak_diff_all/diff_kmeans_no24h/cluster1_no_tcf1.bed',
            quote = F,row.names = F,col.names = F,sep = '\t')
for(i in 1:7){
  hits2 <- findOverlaps(GRanges(CTCF_cluster),GRanges(diff_peak_kmeans[diff_peak_kmeans$cluster_idx==i,]))
  print(length(hits2@from)/hits2@nRnode)
  print(length(hits2@from))
}

CTCF_great <- read.csv('/media/shaoqizhu/easystore/CD8-HP/CTCF/peak_diff_all/diff_kmeans_no24h/cluster1_greatExportAll.tsv',sep = '\t')
great_genes <- unique(unlist(str_split(CTCF_great$Genes[1:20],',')))
for(i in 1:8){
  #print(length(intersect(mm9_deg$gene_name[mm9_deg$cluster==i],great_genes))/length(mm9_deg$gene_name[mm9_deg$cluster==i]))
  print(length(intersect(mm9_deg$gene_name[mm9_deg$cluster==i],great_genes)))
}
library(readxl)
RNA_seq <- read_excel("/media/shaoqizhu/easystore/CD8-HP/RNA_Seq/CD8-HP_CuffDiff_Summary202005.xlsx")
expr <- as.data.frame(RNA_seq[,c(2,3,4,6,5,7,44,45,46,19,20,21)])
rownames(expr) <- RNA_seq$gene_id
expr_great <- expr[rownames(expr) %in% as.character(great_genes),]
plot(log(expr_great$DKO_0h_0),log(expr_great$WT_0h_0))
plot(log(expr_great$WT_72h_0),log(expr_great$WT_0h_0))


PATH='/media/shaoqizhu/easystore/CD8-HP/CTCF/bw_raw/'
name <- c("DKO_24h_1","DKO_24h_2","DKO_72h_1","DKO_72h_2","DKO_na_1","DKO_na_2",
          "WT_24h_1","WT_24h_2","WT_72h_1","WT_72h_2","WT_na_1","WT_na_2")
CTCF_bw <- list()
for (i in 1:12){
  CTCF_bw[[i]] <- import(paste(PATH,name[i],'_treat_pileup.bw',sep = ''))
  CTCF_bw[[i]]$score <- CTCF_bw[[i]]$score/apply(coverage,2,sum)[i]*10^6
  export.bw(CTCF_bw[[i]],paste(PATH,name[i],'_scaled.bw',sep = ''))
}
apply(coverage,2,sum)*10^6

overlap <- matrix(NA,nrow = 6,ncol = 7)
for (i in 1:6){
  for (j in 1:7){
    hits <- findOverlaps(GRanges(merge[unique(diff),][rownames(peak_reorder)[peak_reorder$cluster==i],]),
                         GRanges(diff_peak_kmeans[diff_peak_kmeans$cluster_idx==j,]))
    overlap[i,j] <- length(unique(hits@from)) 
  }
}
colnames(overlap) <- 1:7;rownames(overlap) <- 1:6
pheatmap(overlap,cluster_rows = F,cluster_cols = F,display_numbers = T,number_format = "%.0f",fontsize=20)
table(peak_reorder$cluster)
for (i in 1:6){print(sum(peak_reorder$cluster==i))}
table(diff_peak_kmeans$cluster_idx)

ctcf_diff_peak <- list()
name <- c("DKO_na-WT_na","DKO_72h-WT_72h","WT_72h-WT_na","DKO_72h-DKO_na","DKO_72h-WT_na","DKO_na-WT_72h")
for (i in 1:6){
  ctcf_diff_peak[[i]] <- read.table(paste('/media/shaoqizhu/easystore/CD8-HP/CTCF/peak_diff_all/pair-wise-no_24h/',name[i],'.bed',sep = ''))
}
for (i in 1:6){colnames(ctcf_diff_peak[[i]])[1:3] <- c('chr','start','end')}
for (i in 1:6){
  hits <- findOverlaps(GRanges(ctcf_diff_peak[[i]][ctcf_diff_peak[[i]]$V7<0,]),GRanges(Tcf1))
  print(c(length(unique(hits@from,))))
}

load("/media/shaoqizhu/easystore/Annotation/chromVARmotifs-master/data/mouse_pwms_v1.rda")
library(motifmatchr)
library(BSgenome.Mmusculus.UCSC.mm9)

peak <- diff_peak_kmeans
motif_idx <- matchMotifs(mouse_pwms_v1, GRanges(peak),
                         genome = BSgenome.Mmusculus.UCSC.mm9)

match_num <- apply(motif_idx@assays@data$motifMatches,2,sum)

match_num[which(as.character(motif_idx@colData[,1])=="Ctcf")]

which(as.character(motif_idx@colData$name)=="Ctcf")

peak_motif <- peak[motif_idx@assays@data$motifMatches[,169],]

for (i in 1:7){
  hits <- findOverlaps(GRanges(merge),GRanges(peak_motif[peak_motif$cluster_idx==i,]))
  print(length(unique(hits@from)))
}
overlap <- matrix(NA,nrow = 6,ncol = 7)
for (i in 1:6){
  for (j in 1:7){
    hits <- findOverlaps(GRanges(merge[unique(diff),][rownames(peak_reorder)[peak_reorder$cluster==i],]),
                         GRanges(diff_peak_kmeans[diff_peak_kmeans$cluster_idx==j,]))
    overlap[i,j] <- length(unique(hits@from)) 
  }
}
colnames(overlap) <- 1:7;rownames(overlap) <- 1:6
library(pheatmap)
norm <- apply(overlap,1,sum)%*%t(apply(overlap,2,sum))/sum(overlap)
pheatmap(overlap/norm,cluster_rows = F,cluster_cols = F,display_numbers = T,number_format = "%.0f",fontsize=20)
table(peak_motif$cluster_idx)



