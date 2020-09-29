library(rtracklayer)
library(readr)
library(stringr)
library(pheatmap)
library(grid)
library(openxlsx)
library(limma)
library(dendextend)

cell_cycle <- remm9_degad_table2('/media/shaoqizhu/easystore/CD8-HP/bedFile/GO_term_cell_cycle.txt')
RNA_seq <- read.xlsx("/media/shaoqizhu/easystore/CD8-HP//RNA_Seq/CD8-HP_CuffDiff_Summary201806.xlsx")
setdiff(RNA_seq$gene_id,mm9$gene_name)

non_sample <- which(names(RNA_seq)%in%c("||_","||"))
idx <- as.numeric()
for (i in 1:6){idx = c(idx,(non_sample[i]-7):non_sample[i])}

unique(cell_cycle_deg)

{n <-which(names(RNA_seq)=='End_1')
WT_na_more_than_dKO_na <- RNA_seq$gene_id[as.numeric(RNA_seq[,n-3])<(-1) & as.numeric(RNA_seq[,n-1])<0.05 & as.numeric(RNA_seq[,n-4])>=1]
WT_na_less_than_dKO_na <- RNA_seq$gene_id[as.numeric(RNA_seq[,n-3])>1 & as.numeric(RNA_seq[,n-1])<0.05 & as.numeric(RNA_seq[,n-5])>=1]
n <-which(names(RNA_seq)=='End_5')
WT_na_more_than_WT_s <- RNA_seq$gene_id[as.numeric(RNA_seq[,n-3])<(-1) & as.numeric(RNA_seq[,n-1])<0.05 & as.numeric(RNA_seq[,n-4])>=1]
WT_na_less_than_WT_s <- RNA_seq$gene_id[as.numeric(RNA_seq[,n-3])>1 & as.numeric(RNA_seq[,n-1])<0.05 & as.numeric(RNA_seq[,n-5])>=1]
n <-which(names(RNA_seq)=='End_2')
dKO_s_more_than_dKO_na <- RNA_seq$gene_id[as.numeric(RNA_seq[,n-3])>1 & as.numeric(RNA_seq[,n-1])<0.05 & as.numeric(RNA_seq[,n-5])>=1]
dKO_s_less_than_dKO_na <- RNA_seq$gene_id[as.numeric(RNA_seq[,n-3])<(-1) & as.numeric(RNA_seq[,n-1])<0.05 & as.numeric(RNA_seq[,n-4])>=1]
n <-which(names(RNA_seq)=='End_4')
dKO_s_more_than_WT_s <- RNA_seq$gene_id[as.numeric(RNA_seq[,n-3])>1 & as.numeric(RNA_seq[,n-1])<0.05 & as.numeric(RNA_seq[,n-5])>=1]
dKO_s_less_than_WT_s <- RNA_seq$gene_id[as.numeric(RNA_seq[,n-3])<(-1) & as.numeric(RNA_seq[,n-1])<0.05 & as.numeric(RNA_seq[,n-4])>=1]
n <-which(names(RNA_seq)=='End_3')
dKO_na_more_than_WT_s <- RNA_seq$gene_id[as.numeric(RNA_seq[,n-3])>1 & as.numeric(RNA_seq[,n-1])<0.05 & as.numeric(RNA_seq[,n-5])>=1]
dKO_na_less_than_WT_s <- RNA_seq$gene_id[as.numeric(RNA_seq[,n-3])<(-1) & as.numeric(RNA_seq[,n-1])<0.05 & as.numeric(RNA_seq[,n-4])>=1]
n <-which(names(RNA_seq)=='End_6')
dKO_s_more_than_WT_na <- RNA_seq$gene_id[as.numeric(RNA_seq[,n-3])>1 & as.numeric(RNA_seq[,n-1])<0.05 & as.numeric(RNA_seq[,n-5])>=1]
dKO_s_less_than_WT_na <- RNA_seq$gene_id[as.numeric(RNA_seq[,n-3])<(-1) & as.numeric(RNA_seq[,n-1])<0.05 & as.numeric(RNA_seq[,n-4])>=1]
}

deg <- unique(c(WT_na_more_than_dKO_na,WT_na_less_than_dKO_na,WT_na_more_than_WT_s,WT_na_less_than_WT_s,
  dKO_s_more_than_dKO_na,dKO_s_less_than_dKO_na,dKO_s_more_than_WT_s,dKO_s_less_than_WT_s,
  dKO_na_more_than_WT_s,dKO_na_less_than_WT_s,dKO_s_more_than_WT_na,dKO_s_less_than_WT_na))

cycle <- intersect(unique(cell_cycle$`Gene/Marker`),gene_negative)
cycle <- intersect(unique(cell_cycle$`Gene/Marker`),gene_positive)

library(readxl)
RNA_seq <- read_excel("/media/shaoqizhu/easystore/CD8-HP/RNA_Seq/CD8-HP_CuffDiff_Summary202005.xlsx")

union_deg <- read_excel("/media/shaoqizhu/easystore/CD8-HP/RNA_Seq/CD8-HP_CuffDiff_Summary202005.xlsx",sheet = 14)
deg <- union_deg$gene_id

RNA_seq_deg <- RNA_seq[c(1,which(RNA_seq$...1%in%deg)),]
colnames(RNA_seq_deg) <- RNA_seq_deg[1,]
RNA_seq_deg <- RNA_seq_deg[-1,]
RNA_seq_deg <- RNA_seq_deg[order(RNA_seq_deg$gene_id),]

gene_reorder$gene_name <- rownames(gene_reorder)
gene_reorder <- gene_reorder[order(gene_reorder$gene_name),]
RNA_seq_deg <- cbind(gene_id= gene_reorder$gene_name,cluster=gene_reorder$cluster,RNA_seq_deg[,-1])
RNA_seq_deg <- RNA_seq_deg[order(RNA_seq_deg$cluster),]
write.xlsx(RNA_seq_deg,'/media/shaoqizhu/easystore/CD8-HP/geneList/pair_wise_deg_cluster.xlsx',colNames = T)

WT_na_lt_DKO_na <- as.data.frame(read_excel("/media/shaoqizhu/easystore/CD8-HP/RNA_Seq/CD8-HP_CuffDiff_Summary202005.xlsx",sheet = 2))$gene_id
WT_na_mt_DKO_na <- as.data.frame(read_excel("/media/shaoqizhu/easystore/CD8-HP/RNA_Seq/CD8-HP_CuffDiff_Summary202005.xlsx",sheet = 3))$gene_id

for (i in 2:13){
  print(nrow(read_excel("/media/shaoqizhu/easystore/CD8-HP/RNA_Seq/CD8-HP_CuffDiff_Summary202005.xlsx",sheet = i)))
}
cbind(1:12,excel_sheets("/media/shaoqizhu/easystore/CD8-HP/RNA_Seq/CD8-HP_CuffDiff_Summary202005.xlsx")[2:13])


sample_name <- c( "WT_0h_5n","WT_0h_6n","WT_0h_7n",
                  "DKO_0h_5n","DKO_0h_6n","DKO_0h_7n",
                  "WT1_CD8_72h","WT2_CD8_72h","Ctrl-3s",
                  "DKO1_72h","DKO2_72h","dKO-3s")
condition <- c('WT_na','WT_na','WT_na','DKO_na','DKO_na','DKO_na','WT_s','WT_s','WT_s','DKO_s','DKO_s','DKO_s')
samples <- data.frame(condition);rownames(samples) <- sample_name;samples$condition <- factor(samples$condition)

mtx <- as.data.frame(RNA_seq[RNA_seq$gene_id %in% deg,sample_name])
row.names(mtx) <- RNA_seq$gene_id[RNA_seq$gene_id %in% deg]

mtx_Zscore <- (mtx-apply(X = mtx, MARGIN = 1, mean))/apply(X = mtx, MARGIN = 1, sd)
mtx_hclust <- hclust(dist(mtx_Zscore), method = "complete")
gene_col <- data.frame(cutree(tree = as.dendrogram(mtx_hclust), k = 9))
names(gene_col)<- 'cluster';gene_col$cluster <- factor(gene_col$cluster)

pheatmap(mtx_Zscore,cluster_rows=T,cluster_cols = F,show_rownames = F,cutree_rows=9,
         angle_col='315',annotation_row = gene_col, annotation_col = samples,main = 'Pair-wise DEGs z-score (hclust)')

#fviz_nbclust(mtx_Zscore, kmeans, method = "wss")

set.seed(2);mtx_kmeans <- kmeans(mtx_Zscore,centers=8,nstart = 25)
gene_kmeans <- as.data.frame(sort(mtx_kmeans$cluster)); names(gene_kmeans) <- 'cluster'; 
gene_kmeans$cluster <- factor(gene_kmeans$cluster)

pheatmap(mtx_Zscore[rownames(gene_kmeans),],cluster_rows=F,cluster_cols = F,show_rownames = F,
         angle_col='315',annotation_row = gene_kmeans, annotation_col = samples,main = 'Pair-wise DEGs z-score (k-means)')

reorder_gene <- c(7,6,3,5,1,2,8,4);gene_reorder <- gene_kmeans
for (i in 1:8) {gene_reorder$cluster[gene_kmeans$cluster==reorder_gene[i]] <- i}
name_ordered <- rownames(gene_kmeans)[order(gene_reorder$cluster)]
gene_reorder <- as.data.frame(sort(gene_reorder$cluster));names(gene_reorder)<- 'cluster'
rownames(gene_reorder) <- name_ordered
gene_reorder$cluster <- factor(gene_reorder$cluster)

cbind(1:8,table(gene_reorder$cluster))

pheatmap(mtx_Zscore[rownames(gene_reorder),],cluster_rows=F,cluster_cols = F,show_rownames = F,
         angle_col='315',annotation_row = gene_reorder, annotation_col = samples,main = 'Pair-wise DEGs z-score (k-means)')

write.csv(gene_reorder,'/media/shaoqizhu/easystore/CD8-HP/geneList/gene_kmeans.csv')
mm9_reduced <- read.table('/media/shaoqizhu/easystore/Annotation/mm9_2014_reduced.bed')
names(mm9_reduced) <- c('chr','start','end','width','gene_name','strand')

gene_reorder$gene_name <- rownames(gene_reorder)
mm9_deg <- mm9_reduced[mm9_reduced$gene_name%in%deg,]
gene_dup <- names(table(mm9_deg$gene_name)[table(mm9_deg$gene_name)==2])
gene_reorder <- rbind(gene_reorder,gene_reorder[gene_dup,])
gene_reorder <- gene_reorder[order(gene_reorder$gene_name),]
mm9_deg$cluster <- gene_reorder$cluster

mm9_deg_50kb <- mm9_deg
mm9_deg_50kb$start <- sapply(mm9_deg$start,function(x) ifelse(x>50000,x-50000,0))
mm9_deg_50kb$end <- sapply(mm9_deg$end,function(x) x+50000)

mm9_reduced_50kb <- mm9_reduced
mm9_reduced_50kb$start <- sapply(mm9_reduced$start,function(x) ifelse(x>50000,x-50000,0))
mm9_reduced_50kb$end <- sapply(mm9_reduced$end,function(x) x+50000)

{deg_cluster3 <- mm9_deg_50kb[mm9_deg_50kb$cluster==4,]
peak_cluster3 <- findOverlaps(GRanges(deg_cluster3),GRanges(merge))

peak <- RPKM[unique(peak_cluster3@to),]
peak_Zscore <- (peak-apply(X = peak, MARGIN = 1, mean))/apply(X = peak, MARGIN = 1, sd)
rename <- c("WT_na1","WT_na2","11_WT_na","12_dKO_na","dKO_na1","5_WT_s","WT_s2","10_dKO_s","8_dKO_s","dKO_s1",
            "WT_na","WT_na","WT_na","dKO_na","dKO_na","dKO_s","WT_s","WT_s","dKO_s","dKO_s")

peak_Zscore <- peak_Zscore[,rename[1:10]]

peak_hclust <- hclust(dist(peak_Zscore), method = "complete")
peak_col <- data.frame(cutree(tree = as.dendrogram(peak_hclust), k = 3))

names(peak_col)<- 'cluster';peak_col$cluster <- factor(peak_col$cluster)

sample_DNase <- data.frame(rename[11:20]);names(sample_DNase) <- 'condition';
rownames(sample_DNase) <- rename[1:10];sample_DNase$condition <- factor(sample_DNase$condition)

pheatmap(peak_Zscore,cluster_rows=T,cluster_cols = F,show_rownames = F,cutree_rows=,
         angle_col='315',annotation_row = peak_col, annotation_col = sample_DNase,main = 'Pair-wise diff peaks z-score (hclust)')

set.seed(2);matrix_kmeans <- kmeans(peak_Zscore,centers=3,nstart=25)
peak_kmeans <- as.data.frame(sort(matrix_kmeans$cluster)); names(peak_kmeans) <- 'cluster'; 
peak_kmeans$cluster <- factor(peak_kmeans$cluster)

pheatmap(peak_Zscore[names(sort(matrix_kmeans$cluster)),],cluster_rows=F,cluster_cols = F,show_rownames = F,
         angle_col='315',annotation_row = peak_kmeans, annotation_col = sample_DNase,main = 'Cluster4 DEG associated union peaks (k-means)')
}

mm9_deg_prom <- mm9_deg
for (i in 1:nrow(mm9_deg)){
  if(mm9_deg$strand[i]=='+'){
    ifelse(mm9_deg$start[i]>1000,
           mm9_deg_prom$start[i] <- mm9_deg$start[i]-1000,
           mm9_deg_prom$start[i] <- 0)
    mm9_deg_prom$end[i] <- mm9_deg$start[i]+1000
  }
  if(mm9_deg$strand[i]=='-'){
    ifelse(mm9_deg$end[i]>1000,
           mm9_deg_prom$start[i] <- mm9_deg$end[i]-1000,
           mm9_deg_prom$start[i] <- 0)
    mm9_deg_prom$end[i] <- mm9_deg$end[i]+1000
  }
}


write.table(cbind(mm9_deg),'/media/shaoqizhu/easystore/CD8-HP/bedFile/deg_kmeans.bed',quote = F,col.names = F,row.names = F,sep = '\t')

hits <- suppressWarnings(findOverlaps(GRanges(mm9_deg_50kb),GRanges(diff_peak_kmeans)))
peak_overlap <- cbind(mm9_deg[hits@from,],diff_peak_kmeans[hits@to,])
dist_tss <- sapply(1:nrow(peak_overlap), function(x){peak_overlap[x,9]-peak_overlap[x,2]})
peak_overlap$dist_tss <- dist_tss
peak_overlap <- peak_overlap[order(peak_overlap$cluster,peak_overlap$cluster_idx),]
write.xlsx(peak_overlap,'/media/shaoqizhu/easystore/CD8-HP/geneList/peak_overlap_DNase_RNAseq.xlsx',row.names = F)
write.xlsx(peak_overlap[peak_overlap$cluster==1 & peak_overlap$cluster_idx==5,],'/media/shaoqizhu/easystore/CD8-HP/bedFile/peak_overlap_DNase(5)_RNAseq(1).xlsx',row.names = F)

peak_overlap <- as.data.frame(cbind(peak_overlap$cluster,peak_overlap$cluster_idx))
pheatmap(table(peak_overlap)/norm,display_numbers = TRUE,cluster_cols = F,cluster_rows = F,
         angle_col=0,treeheight_row=10,treeheight_col=10,number_format = "%.2f",fontsize_number = 12,fontsize = 12)

norm <- apply(table(peak_overlap),1,sum)%*%t(apply(table(peak_overlap),2,sum))/sum(table(peak_overlap))

table(peak_reorder$cluster)
hits <- suppressWarnings(findOverlaps(GRanges(mm9_reduced_50kb),GRanges(diff_peak_kmeans)))
table(diff_peak_kmeans[hits@to,5])

cor_para <- round(cor(as.vector(t(correlation)),c(t(table(peak_overlap)/norm))),4)
ggplot()+geom_point(aes(x=as.vector(t(correlation)),y=c(t(table(peak_overlap)/norm))))+
  labs(x='correlation',y='peak overlap',title = paste('cor =',cor_para))+
  theme(plot.title=element_text(size=12,vjust=0.5,hjust=0.5))

for (i in 1:8){print(sum(mtx_kmeans$cluster==i));}
for (i in 1:8){print(length(intersect(unique(cell_cycle$`Gene/Marker`),
                         names(mtx_kmeans$cluster)[mtx_kmeans$cluster==i])))}
for (i in 1:9){print(sum(gene_col$cluster==i));}
for (i in 1:9){print(length(intersect(unique(cell_cycle$`Gene/Marker`),
                                      rownames(gene_col)[gene_col$cluster==i])))}

pca <- as.data.frame(prcomp(mtx_Zscore,scale.=T)$x)
pca$cluster <- as.factor(mtx_kmeans$cluster)
ggplot(data = pca,aes(x=pca$PC1, y=pca$PC2)) + 
  geom_point(aes(color = cluster))# + scale_color_brewer(palette="Dark2")
PCA <- prcomp(mtx_Zscore,scale.=T)
(PCA$sdev/sum(PCA$sdev))[1:2]

overlap <- matrix(data=0,nrow = 8,ncol = 9)
for (i in 1:8){
  for (j in 1:9){
    overlap[i,j] <- length(intersect(rownames(gene_kmeans)[gene_kmeans$cluster==i],
              rownames(gene_col)[gene_col$cluster==j]))
  }
}
overlap <- as.data.frame(overlap)
names(overlap) <- 1:9
pheat <- pheatmap(overlap,display_numbers = TRUE,cluster_cols = F,cluster_rows = F,
                  angle_col=0,treeheight_row=10,treeheight_col=10,number_format = "%.0f",fontsize_number = 12,fontsize = 12)
table(gene_kmeans$cluster)
table(gene_col$cluster)
apply(overlap[pheat$tree_row$order,pheat$tree_col$order],1,sum)

gene_list <- data.frame(matrix(data=NA,nrow=398,ncol=9))
for (i in 1:9){gene_list[1:sum(gene_kmeans$cluster==i),i] <- rownames(gene_kmeans)[gene_kmeans$cluster==i]}
write.csv(gene_kmeans,'/media/shaoqizhu/easystore/CD8_HP_stimulated/gene_kmeans.csv')

{File = '8'
  GO <- read.csv(paste('/media/shaoqizhu/easystore/CD8-HP/DAVID/',File,'.txt',sep=''),header = T, sep = "\t")
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

for (i in 1:6){n[2*i-1] <- which(names(RNA_seq)==paste('End_',i,sep=''))-5;
n[2*i] <- which(names(RNA_seq)==paste('End_',i,sep=''))-4;}
m <- as.character(); for (i in 1:8){m[i] = rownames(gene_kmeans)[which(gene_kmeans$cluster==i)[80]]}
validation <- round(RNA_seq[RNA_seq$gene_id%in%m,n[c(1,2,7,4)]],2)
validation$gene <- RNA_seq$gene_id[as.numeric(rownames(validation))]
rownames(validation) <- order(m)
validation <- validation[order(rownames(validation)),]
pheatmap(validation[,1:4],scale = 'row',cluster_rows = F,cluster_cols = F,angle_col = 315)

great_negative <- read.csv('~/Document/CD8-HP/diff_peaks/great_genes/WT_na-more-than-dKO_na.txt',sep = '\t',header = F)
great_positive <- read.csv('~/Document/CD8-HP/diff_peaks/great_genes/WT_na-less-than-dKO_na.txt',sep = '\t',header = F)
RNA_seq <- read.xlsx("~/Document/CD8-HP/CD8-HP_CuffDiff_Summary201806.xlsx")
n <-which(names(RNA_seq)=='End_1')
WT_na_more_than_dKO_na <- RNA_seq$gene_id[as.numeric(RNA_seq[,n-3])<(-1) & as.numeric(RNA_seq[,n-1])<0.05 & as.numeric(RNA_seq[,n-5])>=1]
WT_na_less_than_dKO_na <- RNA_seq$gene_id[as.numeric(RNA_seq[,n-3])>1 & as.numeric(RNA_seq[,n-1])<0.05 & as.numeric(RNA_seq[,n-4])>=1]

gene_negative <- RNA_seq$gene_id[as.numeric(RNA_seq[,n-3])<(-1) & as.numeric(RNA_seq[,n-1])<0.05 & as.numeric(RNA_seq[,n-5])>=1]
gene_positive <- RNA_seq$gene_id[as.numeric(RNA_seq[,n-3])>1 & as.numeric(RNA_seq[,n-1])<0.05 & as.numeric(RNA_seq[,n-4])>=1]

{Aa = length(intersect(great_negative$V1,gene_negative))
Ab = length(intersect(great_negative$V1,gene_positive))
Ba = length(intersect(great_positive$V1,gene_negative))
Bb = length(intersect(great_positive$V1,gene_positive))

A = length(great_negative$V1)
B = length(great_positive$V1)
a = length(gene_negative)
b = length(gene_positive)}

as.data.frame(cbind(c(Aa,Ab,A),c(Ba,Bb,B),c(a,b,NA)))




