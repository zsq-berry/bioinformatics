library(ggplot2)
library(ggrepel)
library(rtracklayer)
library(stringr)
PATH='/media/shaoqizhu/easystore/CD8-HP/CTCF/narrowPeak/'
name <- c("WT-na","WT-24h","WT-72h","DKO-na","DKO-24h","DKO-72h",
          "WT_na_2","WT_24h_2","WT_72h_2","DKO_na_2","DKO_24h_2","DKO_72h_2")
name1 <- c("WT_na_1","WT_24h_1","WT_72h_1","DKO_na_1","DKO_24h_1","DKO_72h_1",
          "WT_na_2","WT_24h_2","WT_72h_2","DKO_na_2","DKO_24h_2","DKO_72h_2")
CTCF_CD8HP <- list()
for (i in 1:12){CTCF_CD8HP[[i]] <- read.table(paste(PATH,name[i],'_peaks.narrowPeak',sep = ''))}
names(CTCF_CD8HP) <- name1
for (i in 1:12){colnames(CTCF_CD8HP[[i]])[1:3] <- c('chr','start','end')}
for (i in 1:12){CTCF_CD8HP[[i]] <- CTCF_CD8HP[[i]][CTCF_CD8HP[[i]]$V7>=4 & CTCF_CD8HP[[i]]$V9>(-log10(0.05)),]}
for (i in 1:12){write.table(CTCF_CD8HP[[i]],paste(PATH,'Q05_FC4/',name1[i],'_peaks.narrowPeak',sep = ''),
                            quote = F,row.names = F,col.names = F,sep = '\t')}



Venn <- function(n){
  S1 = CTCF_CD8HP[[which(names(CTCF_CD8HP)==n[1])]]
  S2 = CTCF_CD8HP[[which(names(CTCF_CD8HP)==n[2])]]
  S3 = CTCF_CD8HP[[which(names(CTCF_CD8HP)==n[3])]]
  hits1 <- findOverlaps(GRanges(S1),GRanges(S2))
  hits2 <- findOverlaps(GRanges(S1),GRanges(S3))
  hits3 <- findOverlaps(GRanges(S2),GRanges(S3))
  hits0 <- findOverlaps(GRanges(S1[unique(hits1@from),]),GRanges(S3))
  A=nrow(S1);B=nrow(S2);D=nrow(S3);
  c <- length(unique(hits1@from))-length(unique(hits0@from))
  e <- length(unique(hits2@from))-length(unique(hits0@from))
  f <- length(unique(hits3@from))-length(unique(hits0@from))
  g <- length(unique(hits0@from))
  a=A-c-e-g;b=B-c-f-g;d=D-e-f-g
  c(paste(a,b,c,d,e,f,g,A,B,D,sep = ','),
    paste(n[1],"=",A,"','",n[2],"=",B,"','",n[3],"=",D,sep = ''))
}
n=c('DKO_na','DKO_24h','DKO_72h')
Venn(n)

names(CTCF_CD8HP$DKO_na)

n=c('WT_72h','DKO_72h')
{s1 = CTCF_CD8HP[[which(names(CTCF_CD8HP)==n[1])]]
s2 = CTCF_CD8HP[[which(names(CTCF_CD8HP)==n[2])]]
hits <- findOverlaps(GRanges(s1),GRanges(s2))
c(paste(hits@nLnode-length(hits@from),hits@nRnode-length(hits@from),length(hits@from),sep = ','),
paste(n[1],"=",hits@nLnode,"','",n[2],"=",hits@nRnode,sep = ''))
}

overlap <- matrix(NA,nrow = 7,ncol = 7)
for(i in 1:7){
  for (j in 1:7){
    hits <- suppressWarnings(findOverlaps(GRanges(CTCF_CD8HP[[i]]),GRanges(diff_peak_kmeans[diff_peak_kmeans$cluster_idx==j,])))
    overlap[i,j] <- length(unique(hits@from))
  }
}
for (i in 1:7){print(nrow(CTCF_CD8HP[[i]]))}

for (i in 1:7){
  hits1 <- findOverlaps(GRanges(CTCF_CD8HP[[i]]),GRanges(tcf1))
  hits2 <- findOverlaps(GRanges(CTCF_CD8HP[[i]]),GRanges(mm9_promoter))
  hits3 <- findOverlaps(GRanges(CTCF_CD8HP[[i]][hits1@from,]),GRanges(mm9_promoter))
  print(c(length(hits2@from),length(hits3@from)))
}

for (i in 1:7){
  hits <- findOverlaps(GRanges(CTCF_CD8HP[[i]]),GRanges(stat5_dko_24h))
  print(length(hits@from))
}

hits <- findOverlaps(GRanges(CTCF_CD8HP$WT_na),GRanges(Tcf1))
hits1 <- findOverlaps(GRanges(CTCF_CD8HP$WT_na[hits@from,]),GRanges(mm9_promoter))
write.table(CTCF_CD8HP$WT_na[hits@from,][hits1@from,],
            '/media/shaoqizhu/easystore/CD8-HP/CTCF/heatmap/CTCF_WT_na_Tcf1_promoter.bed',
            quote = F,row.names = F,col.names = F,sep = '\t')


hits <- findOverlaps(GRanges(CTCF_CD8HP$WT_na),GRanges(CTCF_CD8HP$DKO_na))
hits1 <- findOverlaps(GRanges(CTCF_CD8HP$WT_na[-hits@from,]),GRanges(tcf1))
WT_na_unique <- CTCF_CD8HP$WT_na[-hits@from,][hits1@from,]
WT_DKO_na <- CTCF_CD8HP$WT_na[hits@from,]

CTCF_bw <- list()
CTCF_bw$DKO_na <- import('/media/shaoqizhu/easystore/CD8-HP/CTCF/bw/DKO-na_treat_pileup.bw')
CTCF_bw$WT_na <- import('/media/shaoqizhu/easystore/CD8-HP/CTCF/bw/WT-na_treat_pileup.bw')

hits2 <- findOverlaps(GRanges(WT_na_unique),CTCF_bw$WT_na)
hits3 <- findOverlaps(GRanges(WT_na_unique),CTCF_bw$DKO_na)
boxplot(log(aggregate(CTCF_bw$WT_na[hits2@to,]$score,list(hits2@from),sum)$x)-log(aggregate(CTCF_bw$DKO_na[hits3@to,]$score,list(hits3@from),sum)$x),
        names=c('WT_na','DKO_na'),main='log coverage on WT_na & DKO_na overlap peaks')



hits <- findOverlaps(GRanges(CTCF_CD8HP$WT_72h),GRanges(diff_peak_kmeans[diff_peak_kmeans$cluster_idx==1,]))

write.table(CTCF_CD8HP$WT_72h[hits@from,1:4],'/media/shaoqizhu/easystore/CD8-HP/CTCF/CTCF_WT_72h_in_cluster1.bed',
            quote = F,row.names = F,col.names = F,sep = '\t')

write.table(diff_peak_kmeans[diff_peak_kmeans$cluster_idx==1,1:4],'/media/shaoqizhu/easystore/CD8-HP/CTCF/cluster1.bed',
            quote = F,row.names = F,col.names = F,sep = '\t')

expr <- as.data.frame(t(RNA_seq[which(RNA_seq$gene_id=='Ctcf'),c(2,3,4,5,6,7,44,45,46,19,20,21)]))
expr$lib <- rownames(expr) 
expr$lib <- factor(expr$lib,expr$lib)
names(expr)[1] <- 'Ctcf'
ggplot(expr,aes(x=lib,y=Ctcf))+geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))





n1 <- c('WT_na',"WT_24h","WT_72h")
n2 <- c('DKO_na',"DKO_24h","DKO_72h")
i=3
{s1 = CTCF_CD8HP[[which(names(CTCF_CD8HP)==n1[i])]]
s2 = CTCF_CD8HP[[which(names(CTCF_CD8HP)==n2[i])]]
hits <- findOverlaps(GRanges(s1),GRanges(s2))
hits1 <- findOverlaps(GRanges(s1[-hits@from,]),GRanges(mm9_promoter))
write.table(s1[-hits@from,][hits1@from,],
            '/media/shaoqizhu/easystore/CD8-HP/CTCF/peak_diff/WT_mt_DKO_in_72h_promoter.bed',
            quote = F,row.names = F,col.names = F,sep = '\t')
write.table(s1[-hits@from,][-hits1@from,],
            '/media/shaoqizhu/easystore/CD8-HP/CTCF/peak_diff/WT_mt_DKO_in_72h_enhancer.bed',
            quote = F,row.names = F,col.names = F,sep = '\t')
}
CTCF_CD8HP$DKO_na[-hits@to,]


PATH='/media/shaoqizhu/easystore/CD8-HP/macs2/narrow_peaks/'
name <- c("WT-na1","WT-na2","11-naive-WT","dKO-na1","dKO-na2","12-navie-dKO",
          "WT-s1","WT-s2","1-stim-WT","5-stim-WT","10-stim-dKO","8-stim-dKO",
          "dKO-s1","dKO-s2")
name1 <- c("WT-na","WT-na","WT-na","dKO-na","dKO-na","dKO-na",
          "WT-s","WT-s","WT-s","WT-s","dKO-s","dKO-s","dKO-s","dKO-s")
DNase_CD8HP <- list()
for (i in 1:length(name)){DNase_CD8HP[[i]] <- read.table(paste(PATH,name[i],'.narrowPeak',sep = ''))}
names(DNase_CD8HP) <- gsub('-','_',name)

for (i in 1:length(name)){colnames(DNase_CD8HP[[i]])[1:3] <- c('chr','start','end')}

DNase_CD8HP_pooled <- list()
for (i in 1:4){
  idx <- which(name1==unique(name1)[i])
  DNase_CD8HP_pooled[[i]] <- DNase_CD8HP[[idx[1]]]
  for (j in idx[2:length(idx)]){
    DNase_CD8HP_pooled[[i]] <- rbind(DNase_CD8HP_pooled[[i]],DNase_CD8HP[[j]])
  }
}
for (i in 1:4){
  DNase_CD8HP_pooled[[i]] <- as.data.frame(reduce(GRanges(DNase_CD8HP_pooled[[i]])))
}
names(DNase_CD8HP_pooled)<- gsub('-','_',unique(name1))

DNase_CD8HP_pooled$WT_na

n=c('WT_na','WT_na')
{
  s1 = CTCF_CD8HP[[which(names(CTCF_CD8HP)==n[1])]]
  s2 = DNase_CD8HP_pooled[[which(names(DNase_CD8HP_pooled)==n[2])]]
  s1 = s1[order(s1$V9,decreasing = T)[1:10000],]
  hits <- findOverlaps(GRanges(s1),GRanges(s2))
  c(paste(hits@nLnode-length(hits@from),hits@nRnode-length(hits@from),length(hits@from),sep = ','),
    paste('CTCF ',n[1],"=",hits@nLnode,"','",'DNase ',n[2],"=",hits@nRnode,sep = ''))
}
nrow(s2)
DNase_CD8HP_pooled$WT_na

which(names(CTCF_CD8HP)==n[1])
DNase_CD8HP$dKO_na1




coverage <- read.csv('/media/shaoqizhu/easystore/CD8-HP/CTCF/merged_peak/no24h/coverage.bed',sep='\t',header = F)
merge <- read.csv('/media/shaoqizhu/easystore/CD8-HP/CTCF/merged_peak/no24h/merged.bed',sep='\t',header = F)
name <- c("DKO_72h","DKO_na","WT_72h","WT_na")
names(coverage)[-1] <- as.character(name)
names(merge)<-c('chr','start','end','idx')
coverage <- coverage[,c(1,5,3,4,2)]

RPK <- coverage[,-1]/coverage[,1]*10^3
RPKM <- as.data.frame(t(t(RPK)/apply(RPK, 2, sum)*10^6))
coverage <- coverage[,-1]

pca <- as.data.frame(prcomp(t(RPKM[,]),scale.=T)$x)
pca$name <- colnames(coverage);pca$exp <- colnames(coverage)

ggplot(data = pca,aes(x=PC1, y=PC2, label=name, colour = exp)) + 
  geom_point() + geom_text_repel(aes(label=name), size = 4, hjust=0, vjust=-1, angle=0)

PCA <- prcomp(t(RPKM[,]),scale.=T)
(PCA$sdev/sum(PCA$sdev))[1:2]

RPKM_FC <- as.data.frame(matrix(NA,nrow = nrow(RPKM),ncol = length(RPKM)*(length(RPKM)-1)/2)); n=0
for (i in 1:(length(RPKM)-1)){
  for (j in (i+1):length(RPKM)){
    n=n+1; RPKM_FC[,n] <- log2(RPKM[,j]/RPKM[,i])
    RPKM_FC[which(RPKM[,j]==0 | RPKM[,i]==0),n] <- 0
    colnames(RPKM_FC)[n] <- paste(colnames(RPKM)[i],colnames(RPKM)[j],sep = '_')
  }
}

diff_idx <- as.numeric()
for (i in 1:length(RPKM_FC)){
  print(c(length(which(RPKM_FC[,i]>=2)),length(which(RPKM_FC[,i]<=-2))))
  diff_idx <- c(diff_idx,which(abs(RPKM_FC[,i])>=2))
  #print(names(RPKM_FC)[i])
}
diff_idx <- unique(diff_idx)
write.table(merge,
            paste('/media/shaoqizhu/easystore/CD8-HP/CTCF/peak_diff/diff_kmeans/merged_peaks.bed',sep = ''),
            quote = F,row.names = F,col.names = F,sep = '\t')

peak_diff <- RPKM[diff_idx,]
peak_diff_z <- (peak_diff-apply(X = peak_diff, MARGIN = 1, mean))/apply(X = peak_diff, MARGIN = 1, sd)

set.seed(2);peak_diff_kmeans <- kmeans(peak_diff_z,centers=5,nstart = 25)
peak_kmeans <- as.data.frame(sort(peak_diff_kmeans$cluster)); names(peak_kmeans) <- 'cluster'; 
peak_kmeans$cluster <- factor(peak_kmeans$cluster)
pheatmap::pheatmap(peak_diff_z[rownames(peak_kmeans),],cluster_rows=F,cluster_cols = F,show_rownames = F,
                   angle_col='315',annotation_row = peak_kmeans, main = 'RNA-seq z-score (k-means)')

reorder_peak <- c(3,2,5,1,4);peak_reorder <- peak_kmeans;
for (i in 1:5) {peak_reorder$cluster[peak_kmeans$cluster==reorder_peak[i]] <- i}
name_ordered <- rownames(peak_kmeans)[order(peak_reorder$cluster)]
peak_reorder <- as.data.frame(sort(peak_reorder$cluster));names(peak_reorder)<- 'cluster'
rownames(peak_reorder) <- name_ordered
peak_reorder$cluster <- factor(peak_reorder$cluster)
pheatmap::pheatmap(peak_diff_z[rownames(peak_reorder),],cluster_rows=F,cluster_cols = F,show_rownames = F,
         angle_col='315',annotation_row = peak_reorder,main = 'CTCF diff peak z-score (k-means)')

tcf1_bw <- import('/media/shaoqizhu/easystore/CD8-HP/CTCF/bw_raw/tcf1.bw')
merge_diff <- merge[diff_idx,][rownames(peak_reorder),]
hits <- findOverlaps(GRanges(merge_diff),tcf1_bw)

tcf1_on_CTCF <- aggregate(tcf1_bw$score[hits@to],list(hits@from),sum)

pheatmap(log(tcf1_on_CTCF$x+0.01),cluster_rows=F,cluster_cols = F,show_rownames = F)

for (i in 1:5){
write.table(merge[rownames(peak_reorder)[peak_reorder$cluster==i],],
            paste('/media/shaoqizhu/easystore/CD8-HP/CTCF/peak_diff/diff_kmeans_no24h/peak_diff_cluster_',i,'.bed',sep = ''),
            quote = F,row.names = F,col.names = F,sep = '\t')
}

Tcf1 <- read.table('/media/shaoqizhu/easystore/CD8-HP/CD8_HP_naive/NaiveCD8_TCF1_BAMPE_peaks_filtered.bed')
colnames(Tcf1)[1:3] <- c('chr','start','end')

for (i in 1:5){
  hits <- findOverlaps(GRanges(Tcf1),GRanges(merge[rownames(peak_reorder)[peak_reorder$cluster==i],]))
  #print(c(length(unique(hits@from)),hits@nRnode))
  print(length(unique(hits@from))/hits@nRnode)
}
CTCF_cluster <- merge[rownames(peak_reorder)[peak_reorder$cluster==4],]
hits <- findOverlaps(GRanges(Tcf1),GRanges(CTCF_cluster))
write.table(CTCF_cluster[-unique(hits@to),],
            '/media/shaoqizhu/easystore/CD8-HP/CTCF/peak_diff/diff_kmeans_no24h/peak_diff_cluster1_no_tcf1.bed',
            quote = F,row.names = F,col.names = F,sep = '\t')
CTCF_great <- read.csv('/media/shaoqizhu/easystore/CD8-HP/CTCF/peak_diff/diff_kmeans_no24h/peak_diff_cluster3_tcf1_greatExportAll.tsv',sep = '\t')
great_genes <- unique(unlist(str_split(CTCF_great$Genes[1:20],',')))
for(i in 1:8){
  #print(length(intersect(mm9_deg$gene_name[mm9_deg$cluster==i],great_genes))/length(mm9_deg$gene_name[mm9_deg$cluster==i]))
  print(intersect(mm9_deg$gene_name[mm9_deg$cluster==i],great_genes))
}

intersect(mm9_deg$gene_name[mm9_deg$cluster==1],great_genes)


hits <- findOverlaps(GRanges(Tcf1),GRanges(CTCF_cluster))
hits1 <- findOverlaps(GRanges(CTCF_cluster[hits@to,]),
                      GRanges(diff_peak_kmeans[diff_peak_kmeans$cluster_idx==1,]))
write.table(CTCF_cluster[hits@to,][hits1@from,],
            '/media/shaoqizhu/easystore/CD8-HP/CTCF/peak_diff/diff_kmeans_no24h/peak_diff_cluster4_in_tcf1_DNase_c1.bed',
            quote = F,row.names = F,col.names = F,sep = '\t')

for(i in 1:7){
  hits2 <- findOverlaps(GRanges(CTCF_cluster),GRanges(diff_peak_kmeans[diff_peak_kmeans$cluster_idx==i,]))
  print(length(hits2@from))
}
library(motifmatchr)
library(BSgenome.Mmusculus.UCSC.mm9)
load("/media/shaoqizhu/easystore/Annotation/chromVARmotifs-master/data/mouse_pwms_v1.rda")
merge <- read.csv('/media/shaoqizhu/easystore/CD8-HP/CTCF/union_peak/merged.bed',sep='\t',header = F)
names(merge)<-c('chr','start','end','idx')
motif_idx <- matchMotifs(mouse_pwms_v1, GRanges(merge),
                         genome = BSgenome.Mmusculus.UCSC.mm9)

match_num <- apply(motif_idx@assays@data$motifMatches,2,sum)

as.data.frame(cbind(as.character(motif_idx@colData[order(match_num,decreasing = T)[1:100],]),
                    sort(match_num,decreasing = T)[1:100]))
match_num[which(as.character(motif_idx@colData[,1])=="Ctcf")]

which(as.character(motif_idx@colData$name)=="Ctcf")

merge_ctcf_motif <- merge[motif_idx@assays@data$motifMatches[,169],]
write.table(merge_ctcf_motif,
            '/media/shaoqizhu/easystore/CD8-HP/CTCF/union_peak/merged_ctcf_motif.bed',
            quote = F,row.names = F,col.names = F,sep = '\t')


for (i in c(3,8,12)){
  print(c(length(which(RPKM_FC[,i]>=2)),length(which(RPKM_FC[,i]<=-2))))
  print(names(RPKM_FC)[i])
  write.table(merge[which(RPKM_FC[,i]>=2),],
              paste('/media/shaoqizhu/easystore/CD8-HP/CTCF/peak_diff/piarwise-diff/',names(RPKM_FC)[i],'_mt_2.bed',sep = ''),
              quote = F,row.names = F,col.names = F,sep = '\t')
  write.table(merge[which(RPKM_FC[,i]<=2),],
              paste('/media/shaoqizhu/easystore/CD8-HP/CTCF/peak_diff/piarwise-diff/',names(RPKM_FC)[i],'_lt_2.bed',sep = ''),
              quote = F,row.names = F,col.names = F,sep = '\t')
}

Tcf1 <- read.table('/media/shaoqizhu/easystore/CD8-HP/CD8_HP_naive/NaiveCD8_TCF1_BAMPE_peaks_filtered.bed')
colnames(Tcf1)[1:3] <- c('chr','start','end')

peak <- merge[which(RPKM_FC$DKO_na_WT_na>=2),]
hits <- findOverlaps(GRanges(Tcf1),GRanges(peak))
hits1 <- findOverlaps(GRanges(peak[-hits@to,]),GRanges(mm9_promoter))
WT_mt_DKO_no_Tcf1_na <- peak[-hits@to,][-hits1@from,]

write.table(tcf1_summit[hits@from,][-hits1@from,],
            '/media/shaoqizhu/easystore/CD8-HP/CTCF/peak_diff/fold_change/WT_mt_DKO_in_Tcf1_summit_na.bed',
            quote = F,row.names = F,col.names = F,sep = '\t')


write.table(diff_peak_kmeans[diff_peak_kmeans$cluster_idx==1,],
            '/media/shaoqizhu/easystore/CD8-HP/CTCF/diff_peak_cluster1.bed',
            quote = F,row.names = F,col.names = F,sep = '\t')


PATH='/media/shaoqizhu/easystore/CD8-HP/CTCF/bw_raw/'
name1<- c("DKO-24h","DKO-72h","DKO-na","WT-24h","WT-72h","WT-na")
CTCF_bw <- list()
for (i in 1:6){
  CTCF_bw[[i]] <- import(paste(PATH,name1[i],'_treat_pileup.bw',sep = ''))
  CTCF_bw[[i]]$score <- CTCF_bw[[i]]$score/apply(coverage,2,sum)[i]*10^6
  export.bw(CTCF_bw[[i]],paste(PATH,name1[i],'_scaled.bw',sep = ''))
}
names(CTCF_bw) <- gsub('-','_',name1)

tcf1_bw <- import(paste(PATH,'tcf1.bw',sep = ''))
tcf1_bw$score <- tcf1_bw$score/5
export.bw(tcf1_bw,paste(PATH,'tcf1_scaled.bw',sep = ''))




CTCF_cluster1_mat <- read.table('/media/shaoqizhu/easystore/CD8-HP/CTCF/heatmap/cluster1_tcf1_regions.bed',header = F)
CTCF_cluster1_region2 <- CTCF_cluster1_mat[CTCF_cluster1_mat$V13=='cluster_1',1:3]
CTCF_cluster1_region2$V4 <- 1:nrow(CTCF_cluster1_region2)
write.table(CTCF_cluster1_region2,
      '/media/shaoqizhu/easystore/CD8-HP/CTCF/heatmap/cluster1_tcf1_regions_c1.bed',
      quote = F,row.names = F,col.names = F,sep = '\t')





library(readxl)
RNA_seq <- read_excel("/media/shaoqizhu/easystore/CD8-HP/RNA_Seq/CD8-HP_CuffDiff_Summary202005.xlsx")
expr <- as.data.frame(RNA_seq[,c(2,3,4,6,5,7,44,45,46,19,20,21)])
rownames(expr) <- RNA_seq$gene_id
CTCF_great <- read.csv('/media/shaoqizhu/easystore/CD8-HP/CTCF/peak_diff/diff_kmeans_no24h/peak_diff_cluster3_tcf1_greatExportAll.tsv',sep = '\t')
great_genes <- unique(unlist(str_split(CTCF_great$Genes[1:20],',')))
#great_genes <- read.csv('/media/shaoqizhu/easystore/CD8-HP/CTCF/peak_diff/WT_mt_DKO_no_Tcf1_na_great.txt',sep = '\t')

expr_great <- expr[rownames(expr) %in% as.character(great_genes),]
expr_z <- (expr_great-apply(X = expr_great, MARGIN = 1, mean))/apply(X = expr_great, MARGIN = 1, sd)
expr_z <- expr_z[-which(expr_z$WT_0h_0=="NaN"),]
pheatmap::pheatmap(expr_z,show_rownames = T,cluster_cols = F,cluster_rows = T)

set.seed(2);expr_kmeans <- kmeans(expr_z,centers=6,nstart = 25)
gene_kmeans <- as.data.frame(sort(expr_kmeans$cluster)); names(gene_kmeans) <- 'cluster'; 
gene_kmeans$cluster <- factor(gene_kmeans$cluster)
pheatmap::pheatmap(expr_z[rownames(gene_kmeans),],cluster_rows=F,cluster_cols = F,show_rownames = F,
                   angle_col='315',annotation_row = gene_kmeans, main = 'RNA-seq z-score (k-means)')

write.csv(rownames(gene_kmeans)[gene_kmeans$cluster==4],
          '/media/shaoqizhu/easystore/CD8-HP/CTCF/peak_diff/geng.csv')

{File = '4'
  GO <- read.csv(paste('/media/shaoqizhu/easystore/CD8-HP/CTCF/peak_diff/david/great_cluster_no_Tcf1_',File,'.txt',sep=''),header = T, sep = "\t")
  GO$Term <- str_split_fixed(GO$Term,'~',2)[,2]
  GO$Term <- factor(GO$Term, levels=rev(unique(GO$Term)))
  GO <-GO[1:20,]
  ggplot(data = GO, aes(x = Count,y = Term),) +
    geom_point(size = 3,aes(color = -log10(FDR))) +
    scale_color_gradient(low="green",high="red") +
    labs(color = expression(-log[10](FDR)),x = "Gene Count",y="GO term",
         title=paste('DAVID GO in cluster',File)) +
    theme(plot.title=element_text(size=12,vjust=0.5,hjust=0.5),
          legend.title=element_text(size=10),
          legend.text=element_text(size=10),
          legend.position="right",
          axis.text=element_text(size=10),
          axis.title=element_text(size=10),)
}

library(GenomicRanges)
findOverlaps(GRanges(CTCF_CD8HP$WT_na_1),GRanges(Tcf1))
