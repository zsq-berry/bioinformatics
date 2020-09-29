#======================================= HiC AB compartment

AB_H3K27ac <- read.table('/media/shaoqizhu/easystore/HiC_CD8HP/eigen_100kb/AB_eigen/naive-WT-CD8-1_201903.bed')
AB_domain <- read.table('/media/shaoqizhu/easystore/HiC_CD8HP/eigen_100kb/AB_eigen/AB_eigen.txt')
colnames(AB_domain)[4:11] <- c("naive_DKO_2016","naive_DKO_2019","naive_WT_2016","naive_WT_2019","stim_DKO_2016","stim_DKO_2019","stim_WT_2016","stim_WT_2019")
AB_domain$V1 <- as.character(AB_domain$V1)
for (i in unique(AB_domain$V1)){
  L <- nrow(AB_domain[AB_domain$V1 %in% i,])
  pm <- (apply(AB_domain[AB_domain$V1 %in% i,][,4:11],2,function(x){cor(x,AB_H3K27ac$V15[AB_H3K27ac$V1 %in% i])}))
  for (j in 1:8){
    if(pm[j]<0)
      AB_domain[AB_domain$V1 %in% i,j+3] <- -AB_domain[AB_domain$V1 %in% i,j+3]
  }
}
cor(AB_domain[AB_domain$V1 %in% 'chr13',][,4:11],AB_H3K27ac$V15[AB_H3K27ac$V1 %in% 'chr13'])
nrow(AB_domain[AB_domain$V1 %in% 'chr1',])

library(pheatmap)
library(RColorBrewer)
library(limma)
library(tidyr)
library(ggplot2)

batch <- c(1,2,1,2,1,2,1,2)
AB_batch <- AB_domain
for (i in unique(AB_domain$V1)){
  AB_batch[AB_batch$V1 %in% i,4:11] <- removeBatchEffect(AB_domain[AB_domain$V1 %in% i,4:11], batch)
  }
AB_compare <- AB_domain[AB_domain$V1 %in% 'chr19',4:11]
pheatmap(cor(AB_batch[,4:11]),cluster_rows = T,cluster_cols = T,
         display_numbers=T,fontsize_number=10,angle_col = 315,number_format = "%.4f",
         colorRampPalette(brewer.pal(n = 7, name = "OrRd"))(100))

plotDat <- gather(AB_compare,'Sample','Eigenvector')
ggplot(plotDat, aes(Sample, Eigenvector)) + geom_boxplot() +labs(title='Batch effect removed') + 
  theme(plot.title=element_text(size=12,vjust=0.5,hjust=0.5),
        axis.text.x = element_text(angle = 315, vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10))

cor(AB_domain$naive_WT_2019,AB_batch$naive_WT_2019)


#====================================== HiC loops

loops_WT_DKO <- read.table('/media/shaoqizhu/easystore/HiC_CD8HP/hiccups/new/loops_WT-DKO_in_naive.bed')
loops_WT_DKO <- unique(loops_WT_DKO)
write.table(cbind(paste(loops_WT_DKO$V1,loops_WT_DKO$V2,loops_WT_DKO$V3,sep = ':'),
                  paste(loops_WT_DKO$V4,loops_WT_DKO$V5,loops_WT_DKO$V6,sep = ':'))
            ,'/media/shaoqizhu/easystore/HiC_CD8HP/hiccups/new/loops_WT-DKO_in_naive_unique.bed',
            quote = F,col.names = F,row.names = F,sep = '\t')

loops_naive_stim <- read.table('/media/shaoqizhu/easystore/HiC_CD8HP/hiccups/new/loops_naive-stim_in_WT.bed')
loops_naive_stim <- unique(loops_naive_stim)
write.table(cbind(paste(loops_naive_stim$V1,loops_naive_stim$V2,loops_naive_stim$V3,sep = ':'),
                  paste(loops_naive_stim$V4,loops_naive_stim$V5,loops_naive_stim$V6,sep = ':'))
            ,'/media/shaoqizhu/easystore/HiC_CD8HP/hiccups/new/loops_naive-stim_in_WT_unique.bed',
            quote = F,col.names = F,row.names = F,sep = '\t')



merged_loop <- read.table('/media/shaoqizhu/easystore/HiC_CD8HP/merged_loop.bed')

merged_loop <- unique(merged_loop)

write.table(cbind(paste(merged_loop$V1,merged_loop$V2,merged_loop$V2,sep = ':'),
                  paste(merged_loop$V1,merged_loop$V3,merged_loop$V3,sep = ':'))
            ,'/media/shaoqizhu/easystore/HiC_CD8HP/merged_loop.bed',
            quote = F,col.names = F,row.names = F,sep = '\t')

cbind(paste(merged_loop$V1,merged_loop$V2,merged_loop$V2,sep = ':'),
paste(merged_loop$V1,merged_loop$V3,merged_loop$V3,sep = ':'))

merged_loop <- read.table('/media/shaoqizhu/easystore/HiC_CD8HP/hiccups/new/loops_naive-stim_in_WT.bed')
merged_loop <- unique(merged_loop)
merged_loop_density <- read.table('/media/shaoqizhu/easystore/HiC_CD8HP/hiccups/new/loops_naive-stim_in_WT/loops_naive-stim_in_WT.txt')
colnames(merged_loop_density) <- c("naive_WT_2016", "naive_WT_2019", "stim_WT_2016", "stim_WT_2019")

merged_loop <- read.table('/media/shaoqizhu/easystore/HiC_CD8HP/hiccups/new/loops_WT-DKO_in_naive.bed')
merged_loop <- unique(merged_loop)
merged_loop_density <- read.table('/media/shaoqizhu/easystore/HiC_CD8HP/hiccups/new/loops_WT-DKO_in_naive/loops_WT-DKO_in_naive.txt')
colnames(merged_loop_density) <- c("naive_DKO_2016", "naive_DKO_2019", "naive_WT_2016", "naive_WT_2019")

RPKM <- as.data.frame(t(t(merged_loop_density)/apply(merged_loop_density, 2, sum)*10^6))

pca <- as.data.frame(prcomp(t(RPKM),scale.=T)$x)
pca$name <- rownames(pca);pca$exp <- paste(str_split_fixed(rownames(pca),'_',3)[,1],str_split_fixed(rownames(pca),'_',3)[,2],sep = '_')


ggplot(data = pca,aes(x=PC1, y=PC2, label=name, colour = exp)) + 
  geom_point() + geom_text_repel(aes(label=name), size = 4, hjust=0, vjust=-1, angle=0)


PCA <- prcomp(t(coverage[,]),scale.=T)
PCA <- prcomp(t(RPKM[,]),scale.=T)
(PCA$sdev/sum(PCA$sdev))[1:2]

diff <- numeric(); N = length(unique(pca$exp))
for (i in 1:(N-1)){
  for (j in (i+1):N){
    condition1 = which(pca$exp==unique(pca$exp)[i])
    condition2 = which(pca$exp==unique(pca$exp)[j])
    group <- factor(c(rep(i,length(condition1)),rep(j,length(condition2))))
    y = DGEList(counts = merged_loop_density[,c(condition1,condition2)],group = group)
    print(names(merged_loop_density[,c(condition1,condition2)]))
    y <- calcNormFactors(y)
    y <- estimateCommonDisp(y)
    fit <- glmQLFit(y)
    qlf <- glmQLFTest(fit)
    fdr <- p.adjust(p = qlf$table$PValue,method = "BH")
    diff_idx <- which(qlf$table$PValue < 0.1 & abs(qlf$table$logFC) > log2(1.5))
    diff <- c(diff,diff_idx)
    change <- cbind(merged_loop[diff_idx,],qlf$table$PValue[diff_idx],fdr[diff_idx],qlf$table$logFC[diff_idx])
    write.table(change,paste('/media/shaoqizhu/easystore/HiC_CD8HP/hiccups/new/loops_WT-DKO_in_naive/',unique(pca$exp)[i],'-',unique(pca$exp)[j],'.bed',sep=''),quote = F,sep = '\t',row.names = F,col.names = F)
    #write.table(change,paste('/media/shaoqizhu/easystore/HiC_CD8HP/hiccups/new/loops_naive-stim_in_WT/',unique(pca$exp)[i],'-',unique(pca$exp)[j],'.bed',sep=''),quote = F,sep = '\t',row.names = F,col.names = F)
  }
}
unique(diff)

#============================= DKO and WT in naive
loop_naive <- read.table('/media/shaoqizhu/easystore/HiC_CD8HP/hiccups/new/loops_WT-DKO_in_naive/naive_DKO-naive_WT.bed')
loop_naive_WT <- loop_naive[loop_naive$V9>0,]
loop_naive_DKO <- loop_naive[loop_naive$V9<0,]

colnames(loop_naive_WT)[1:3] <- c('chr','start','end'); colnames(loop_naive_WT)[4:6] <- c('chr','start','end');
colnames(loop_naive_DKO)[1:3] <- c('chr','start','end'); colnames(loop_naive_DKO)[4:6] <- c('chr','start','end');


anchor_naive_WT <- rbind(loop_naive_WT[,1:3],loop_naive_WT[,4:6])
anchor_naive_DKO <- rbind(loop_naive_DKO[,1:3],loop_naive_DKO[,4:6])
anchor_naive_WT$chr <- paste('chr',anchor_naive_WT$chr,sep = '')
anchor_naive_DKO$chr <- paste('chr',anchor_naive_DKO$chr,sep = '')

WT_mt_DKO <- read_excel("/media/shaoqizhu/easystore/CD8-HP/RNA_Seq/CD8-HP_CuffDiff_Summary202005.xlsx",sheet = 3)
WT_lt_DKO <- read_excel("/media/shaoqizhu/easystore/CD8-HP/RNA_Seq/CD8-HP_CuffDiff_Summary202005.xlsx",sheet = 2)

hits1 <- findOverlaps(GRanges(anchor_naive_WT),GRanges(mm9_promoter[mm9_promoter$gene_name %in% WT_mt_DKO$gene_id,]))
hits2 <- findOverlaps(GRanges(anchor_naive_DKO),GRanges(mm9_promoter[mm9_promoter$gene_name %in% WT_mt_DKO$gene_id,]))
hits3 <- findOverlaps(GRanges(anchor_naive_WT),GRanges(mm9_promoter[mm9_promoter$gene_name %in% WT_lt_DKO$gene_id,]))
hits4 <- findOverlaps(GRanges(anchor_naive_DKO),GRanges(mm9_promoter[mm9_promoter$gene_name %in% WT_lt_DKO$gene_id,]))
cbind(c(length(unique(hits1@to)),length(unique(hits2@to))),
      c(length(unique(hits3@to)),length(unique(hits4@to))))

mm9_promoter[mm9_promoter$gene_name %in% WT_mt_DKO$gene_id,][hits1@to,5]

colnames(merged_loop)[1:3] <- c('chr','start','end'); colnames(merged_loop)[4:6] <- c('chr','start','end');
anchor_merged_loop <- rbind(merged_loop[,1:3],merged_loop[,4:6])
anchor_merged_loop$chr <- paste('chr',anchor_merged_loop$chr,sep = '')
hits <- findOverlaps(GRanges(anchor_merged_loop),
                     GRanges(mm9_promoter[mm9_promoter$gene_name %in% WT_mt_DKO$gene_id,]))

gene_loops <- mm9_promoter[mm9_promoter$gene_name %in% WT_mt_DKO$gene_id,][hits@to,5]

idx = hits@from[which(gene_loops=='Myb')]
qlf$table$PValue[idx];fdr[idx];qlf$table$logFC[idx]


length((hits2@from))
hits <- hits3;L=nrow(loop_naive_WT);anchor <- 0
hits <- hits4;L=nrow(loop_naive_DKO);anchor <- 0
for (i in 1:length(hits@from)){if (hits@from[i] > L){anchor[i] <- hits@from[i]-L}else{anchor[i] <- hits@from[i]}}
length(unique(anchor))

#============================= naive and stim in WT
loop_WT <- read.table('/media/shaoqizhu/easystore/HiC_CD8HP/hiccups/new/loops_naive-stim_in_WT/naive_WT-stim_WT.bed')
loop_WT_naive <- loop_WT[loop_WT$V9<0,]
loop_WT_stim <- loop_WT[loop_WT$V9>0,]

colnames(loop_WT_naive)[1:3] <- c('chr','start','end'); colnames(loop_WT_naive)[4:6] <- c('chr','start','end');
colnames(loop_WT_stim)[1:3] <- c('chr','start','end'); colnames(loop_WT_stim)[4:6] <- c('chr','start','end');


anchor_WT_naive <- rbind(loop_WT_naive[,1:3],loop_WT_naive[,4:6])
anchor_WT_stim <- rbind(loop_WT_stim[,1:3],loop_WT_stim[,4:6])
anchor_WT_naive$chr <- paste('chr',anchor_WT_naive$chr,sep = '')
anchor_WT_stim$chr <- paste('chr',anchor_WT_stim$chr,sep = '')

naive_mt_stim <- read_excel("/media/shaoqizhu/easystore/CD8-HP//RNA_Seq/CD8-HP_CuffDiff_Summary202005.xlsx",sheet = 13)
naive_lt_stim <- read_excel("/media/shaoqizhu/easystore/CD8-HP//RNA_Seq/CD8-HP_CuffDiff_Summary202005.xlsx",sheet = 12)

hits1 <- findOverlaps(GRanges(anchor_WT_naive),GRanges(mm9_promoter[mm9_promoter$gene_name %in% naive_mt_stim$gene_id,]))
hits2 <- findOverlaps(GRanges(anchor_WT_stim),GRanges(mm9_promoter[mm9_promoter$gene_name %in% naive_mt_stim$gene_id,]))
hits3 <- findOverlaps(GRanges(anchor_WT_naive),GRanges(mm9_promoter[mm9_promoter$gene_name %in% naive_lt_stim$gene_id,]))
hits4 <- findOverlaps(GRanges(anchor_WT_stim),GRanges(mm9_promoter[mm9_promoter$gene_name %in% naive_lt_stim$gene_id,]))
cbind(c(length(unique(hits1@to)),length(unique(hits2@to))),
      c(length(unique(hits3@to)),length(unique(hits4@to))))

hits <- hits3;L=nrow(loop_WT_naive);anchor <- 0
hits <- hits4;L=nrow(loop_WT_stim);anchor <- 0
for (i in 1:length(hits@from)){if (hits@from[i] > L){anchor[i] <- hits@from[i]-L}else{anchor[i] <- hits@from[i]}}
length(unique(anchor))

length(unique(hits4@to))
#================================
plot(merged_loop_density$stim_DKO_2016,merged_loop_density$stim_DKO_2019)

for (i in 1:19){
  chr_eigen <- gsub(' ','',cbind(as.character(mm9_size$V1[i]),format(seq(0,mm9_size$V2[i],10^5), scientific=F),
                                 c(format(seq(0,mm9_size$V2[i],10^5), scientific=F),mm9_size$V2[i])[-1]))
  write.table(chr_eigen,paste('/media/shaoqizhu/easystore/HiC_CD8HP/eigen_100kb/chrome/',mm9_size$V1[i],'.txt',sep = ''),
              quote = F,col.names = F,row.names = F,sep = '\t')
}
mm9_size$V1[1]

loop_naive_WT <- read.table('/media/shaoqizhu/easystore/HiC_CD8HP/hiccups/hiccupDiff/naive-stim/differential_loops1.bedpe')
loop_stim_WT <- read.table('/media/shaoqizhu/easystore/HiC_CD8HP/hiccups/hiccupDiff/naive-stim/differential_loops2.bedpe')
colnames(loop_naive_WT)[1:3] <- c('chr','start','end'); colnames(loop_naive_WT)[4:6] <- c('chr','start','end');
colnames(loop_stim_WT)[1:3] <- c('chr','start','end'); colnames(loop_stim_WT)[4:6] <- c('chr','start','end');
which(loop_naive_WT$V18<10^5)


anchor_naive_WT <- rbind(loop_naive_WT[,1:3],loop_naive_WT[,4:6])
anchor_stim_WT <- rbind(loop_stim_WT[,1:3],loop_stim_WT[,4:6])
anchor_naive_WT$chr <- paste('chr',anchor_naive_WT$chr,sep = '')
anchor_stim_WT$chr <- paste('chr',anchor_stim_WT$chr,sep = '')
hits1 <- findOverlaps(GRanges(anchor_naive_WT),GRanges(mm9_deg_prom))
hits2 <- findOverlaps(GRanges(anchor_stim_WT),GRanges(mm9_deg_prom))
cbind(table(mm9_deg_prom[hits1@to,]$cluster),
      table(mm9_deg_prom[hits1@to,]$cluster)/table(mm9_deg_prom$cluster),
      table(mm9_deg_prom[hits2@to,]$cluster),
      table(mm9_deg_prom[hits2@to,]$cluster)/table(mm9_deg_prom$cluster),
      table(mm9_deg_prom$cluster))

hits <- findOverlaps(GRanges(anchor_stim_WT[hits2@from,]),GRanges(ctcf))
length(unique(hits@from))

library(readxl)
naive_mt_stim <- read_excel("/media/shaoqizhu/easystore/CD8-HP//RNA_Seq/CD8-HP_CuffDiff_Summary202005.xlsx",sheet = 13)
naive_lt_stim <- read_excel("/media/shaoqizhu/easystore/CD8-HP//RNA_Seq/CD8-HP_CuffDiff_Summary202005.xlsx",sheet = 12)
mm9_promoter[mm9_promoter$gene_name %in% naive_mt_stim$gene_id,]
hits1 <- findOverlaps(GRanges(anchor_naive_WT),
                      GRanges(mm9_promoter[mm9_promoter$gene_name %in% naive_mt_stim$gene_id,]))
hits2 <- findOverlaps(GRanges(anchor_stim_WT),
                      GRanges(mm9_promoter[mm9_promoter$gene_name %in% naive_mt_stim$gene_id,]))


cbind(anchor_stim_WT[hits2@from,],mm9_promoter[mm9_promoter$gene_name %in% naive_lt_stim$gene_id,][hits2@to,])
loop_stim_WT[340,]

for (i in 1:8){
hits3 <- findOverlaps(GRanges(anchor_stim_WT[hits2@from,][which(mm9_deg_50kb[hits2@to,]$cluster==i),]),GRanges(ctcf))
print(length(unique(hits3@from)))
}




cbind(anchor_naive_DKO[hits2@from,],mm9_promoter[mm9_promoter$gene_name %in% WT_lt_DKO$gene_id,][hits2@to,])
loop_stim_WT[340,]

for (i in 1:8){
  hits3 <- findOverlaps(GRanges(anchor_stim_WT[hits2@from,][which(mm9_deg_50kb[hits2@to,]$cluster==i),]),GRanges(ctcf))
  print(length(unique(hits3@from)))
}

#================================== TAD
#---------------------------------- TAD_naive
merged_tad <- read.table('/media/shaoqizhu/easystore/HiC_CD8HP/TAD_naive/WT_CD8_naive_pooled_TAD.bedpe')
merged_tad <- merged_tad[,1:3]
write.table(cbind(paste(merged_tad$V1,merged_tad$V2,merged_tad$V3,sep = ':'),
                  paste(merged_tad$V1,merged_tad$V2,merged_tad$V3,sep = ':')),
            '/media/shaoqizhu/easystore/HiC_CD8HP/TAD_naive/merged_tad_straw.bed',quote = F,col.names = F,row.names = F,sep = '\t')
merged_tad$V1 <- paste('chr',merged_tad$V1,sep = '')
colnames(merged_tad) <- c('chr','start','end')

#---------------------------------- end
merged_tad <- read.table('/media/shaoqizhu/easystore/HiC_CD8HP/TAD/merged_tad.bed')
write.table(cbind(paste(merged_tad$V1,merged_tad$V2,merged_tad$V3,sep = ':'),
                  paste(merged_tad$V1,merged_tad$V2,merged_tad$V3,sep = ':')),
            '/media/shaoqizhu/easystore/HiC_CD8HP/TAD/merged_tad_straw.bed',quote = F,col.names = F,row.names = F,sep = '\t')
merged_tad$V1 <- paste('chr',merged_tad$V1,sep = '')
colnames(merged_tad) <- c('chr','start','end')

mm9_size <- read.table('/media/shaoqizhu/easystore/read-through/annotation/mm9.genome')

mm9_size$V3 <- floor(mm9_size$V2/10^4)
merged_tad$width <- (merged_tad$end-merged_tad$start)/10^4

set.seed(2)
random_tad <- as.numeric()
random_chr <- rep(merged_tad$chr,each=100)
for (i in 1:nrow(merged_tad)){
  start <- sample(30:(mm9_size$V3[mm9_size$V1 == paste('chr',merged_tad$chr[i],sep = '')]-merged_tad$width[i]),100)
  random_tad <- rbind(random_tad,cbind(start*10^4,(start+merged_tad$width[i])*10^4))
}

random_tad <- gsub(' ','',format(random_tad,scientific=F))

random_tad_straw <- cbind(paste(random_chr,random_tad[,1],random_tad[,2],sep = ':'),
                          paste(random_chr,random_tad[,1],random_tad[,2],sep = ':'))
write.table(random_tad_straw,
            '/media/shaoqizhu/easystore/HiC_CD8HP/TAD/random_tad_straw.bed',quote = F,col.names = F,row.names = F,sep = '\t')

random_tad_pet <- read.table('/media/shaoqizhu/easystore/HiC_CD8HP/TAD/random_tad/random_tad_pet.txt')
random_tad_pet$idx <- rep(1:1472,each=120)
random_tad_mean <- aggregate(random_tad_pet[,1:8],list(random_tad_pet$idx),mean)[,-1]
colnames(random_tad_mean) <- c("naive_DKO_2016","naive_DKO_2019","naive_WT_2016","naive_WT_2019","stim_DKO_2016","stim_DKO_2019","stim_WT_2016","stim_WT_2019")
random_tad_mean <- as.data.frame(t(t(random_tad_mean)/apply(random_tad_mean, 2, sum)*10^7))


#merged_tad_pet <- normalizeQuantiles(merged_tad_pet-sum(merged_tad_pet)/8/1472)
# groups <- factor(c(1,1,2,2,3,3,4,4))
# AOV <- rep(0,1472)
# for(i in 1:1472){AOV[i] <- oneway.test(as.numeric(merged_tad_pet[i,])~groups)$p.value}
# which(AOV<0.5)

merged_tad_pet <- read.table('/media/shaoqizhu/easystore/HiC_CD8HP/TAD/merged_tad/merged_tad_pet.txt')
colnames(merged_tad_pet)[1:8] <- c("naive_DKO_2016","naive_DKO_2019","naive_WT_2016","naive_WT_2019","stim_DKO_2016","stim_DKO_2019","stim_WT_2016","stim_WT_2019")
merged_tad_pet <- as.data.frame(t(t(merged_tad_pet)/apply(merged_tad_pet, 2, sum)*10^7))
merged_tad_pet <- merged_tad_pet/random_tad_mean
merged_tad_pet <- normalizeQuantiles(merged_tad_pet-sum(merged_tad_pet)/8/1472)
#merged_tad_pet <- merged_tad_pet[,c(3,1,7,5)]
#merged_tad_pet <- merged_tad_pet[,c(4,2,8,6)]

merged_tad_pet_Z <- (merged_tad_pet-apply(X = merged_tad_pet, MARGIN = 1, mean))/apply(X = merged_tad_pet, MARGIN = 1, sd)

set.seed(2);TAD_kmeans <- kmeans(merged_tad_pet_Z,centers=5,nstart=25)
TAD_kmeans$reorder <- TAD_kmeans$cluster; reorders <- c(1,2,3,4,5)# 2019 reorders <- c(5,1,4,2,3);  2016 reorders <- c(4,2,3,5,1)
for (i in 1:5){TAD_kmeans$reorder[TAD_kmeans$cluster==i] <- reorders[i]}
tad_kmeans <- as.data.frame(sort(TAD_kmeans$reorder)); names(tad_kmeans) <- 'cluster'; 
tad_kmeans$cluster <- factor(tad_kmeans$cluster)
rownames(tad_kmeans) <- order(TAD_kmeans$reorder)
pheatmap(merged_tad_pet_Z[order(TAD_kmeans$reorder),], cluster_rows = F,cluster_cols = F,show_rownames = F,
         angle_col='315',annotation_row = tad_kmeans)


gene <- mm9_reduced[grep('Cdk1',mm9_reduced$gene_name),]
gene <- mm9_reduced[mm9_reduced$gene_name=='Rb1',]
hits <- findOverlaps(GRanges(merged_tad),GRanges(gene))

merged_tad[hits@from,]
merged_tad_pet[hits@from,]
TAD_kmeans$cluster[hits@from]

merged_tad_pet[TAD_kmeans$cluster==5,][23,]
merged_tad[TAD_kmeans$cluster==5,][23,]

gene_FC <- apply(mtx[,7:9],1,mean)/apply(mtx[,1:3],1,mean)

hits <- findOverlaps(GRanges(merged_tad),GRanges(mm9_deg))
cbind(merged_tad[hits@from,],TAD_cluster=TAD_kmeans$cluster[hits@from],mm9_deg[hits@to,])[mm9_deg[hits@to,]$cluster%in%c(1:3),]


overlap <- matrix(data = NA, nrow = 5, ncol = 7)
for (i in 1:5){
  for (j in 1:7){
    overlap[i,j] <- length(findOverlaps(GRanges(merged_tad[TAD_kmeans$reorder==i,]),
             GRanges(diff_peak_kmeans[diff_peak_kmeans$cluster_idx==j,]))@from)
}
}
rownames(overlap) <- 1:5; colnames(overlap) <- 1:7
norm <- apply(overlap,1,sum)%*%t(apply(overlap,2,sum))/sum(overlap)
#table(TAD_kmeans$reorder)%*%t(table(diff_peak_kmeans$cluster_idx))/10^3
pheatmap(overlap/norm,cluster_cols = F,cluster_rows = F,display_numbers = T,
         angle_col=0,number_format = "%.2f",fontsize = 20)


overlap <- matrix(data = NA, nrow = 5, ncol = 8)
for (i in 1:5){
  for (j in 1:8){
    overlap[i,j] <- length(findOverlaps(GRanges(merged_tad[TAD_kmeans$reorder==i,]),
                                        GRanges(mm9_deg[mm9_deg$cluster==j,]))@from)
  }
}
rownames(overlap) <- 1:5; colnames(overlap) <- 1:8
norm <- apply(overlap,1,sum)%*%t(apply(overlap,2,sum))/sum(overlap)
pheatmap(overlap/norm,cluster_cols = F,cluster_rows = F,display_numbers = T,
         angle_col=0,number_format = "%.2f",fontsize = 20)



#================================== Tcf1 & TAD
#bedtools coverage -a naive_tad.bed -b /media/shaoqizhu/easystore/CD8-HP/CD8_HP_naive/NaiveCD8_TCF1_BAMPE_peaks_filtered.bed > naive_tad_tcf1.bed
merged_tad_pet=as.data.frame(intra_tad/total_tad)
colnames(merged_tad_pet) <- c("naive_DKO_2016","naive_DKO_2019","naive_WT_2016","naive_WT_2019")
merged_tad_tcf1 <- read.table('/media/shaoqizhu/easystore/HiC_CD8HP/TAD_naive/naive_tad_tcf1.bed')
merged_tad_tcf1$V8 = merged_tad_tcf1$V4/merged_tad_tcf1$V6
tad_tcf1_high <- order(merged_tad_tcf1$V8,decreasing = T)[1:200]
tad_tcf1_low <- order(merged_tad_tcf1$V8,decreasing = F)[1:200]
tad_high_ratio <-(merged_tad_pet$naive_DKO_2016[tad_tcf1_high]/merged_tad_pet$naive_WT_2016[tad_tcf1_high])
tad_low_ratio <- (merged_tad_pet$naive_DKO_2016[tad_tcf1_low]/merged_tad_pet$naive_WT_2016[tad_tcf1_low])
boxplot(tad_high_ratio,tad_low_ratio,names=c('high','low'),main='WT TAD (stim/naive)')
boxplot(tad_high_ratio,tad_low_ratio,names=c('Tcf1_high','Tcf1_low'),main='naive 2016 TAD (DKO-WT)')

df <- data.frame(
  x = c(tad_high_ratio,tad_low_ratio),
  Tcf1_coverage = factor(c(rep('high',length(tad_high_ratio)),rep('low',length(tad_low_ratio))))
)

ggplot(df, aes(x,colour=Tcf1_coverage)) + stat_ecdf(geom = "step") + xlab('WT TAD log(stim/naive)')
ggplot(df, aes(x,colour=Tcf1_coverage)) + stat_ecdf(geom = "step") + xlab('naive TAD log(DKO/WT)')
ks.test(tad_high_ratio,tad_low_ratio)


boxplot(merged_tad_pet[tad_tcf1_high,]$naive_WT_2016,
        merged_tad_pet[tad_tcf1_low,]$naive_WT_2016,
        names=c('Tcf1_high','Tcf1_low'),main='TAD score in naive-WT 2016')

merged_tad[tad_tcf1_low[1:10],]
merged_tad_pet[tad_tcf1_low[1:10],c(5,7)]
merged_tad_tcf1[tad_tcf1_high[1:10],]
table(TAD_kmeans$reorder[tad_tcf1_high])

sort(merged_tad_tcf1$V7,decreasing = T)[1:100]
merged_tad_stat5 <- read.table('/media/shaoqizhu/easystore/HiC_CD8HP/TAD/merged_tad_stat5_peak.bed')
tad_stat5_high <- order(merged_tad_stat5$V7,decreasing = T)[1:200]
tad_stat5_low <- order(merged_tad_stat5$V7,decreasing = F)[1:200]
tad_high_ratio <- (merged_tad_pet$naive_DKO_2016[tad_stat5_high]-merged_tad_pet$naive_WT_2016[tad_stat5_high])
tad_low_ratio <- (merged_tad_pet$naive_DKO_2016[tad_stat5_low]-merged_tad_pet$naive_WT_2016[tad_stat5_low])
boxplot(tad_high_ratio,tad_low_ratio,names=c('high','low'),main='WT TAD log(stim/naive)')
boxplot(tad_high_ratio,tad_low_ratio,names=c('high','low'),main='naive TAD log(DKO/WT)')

mean(sort(merged_tad_stat5$V7,decreasing = T)[1:100])/mean(sort(merged_tad_stat5$V7,decreasing = F)[1:100])

df <- data.frame(
  x = c(tad_high_ratio,tad_low_ratio),
  Stat5_coverage = factor(c(rep('high',length(tad_high_ratio)),rep('low',length(tad_low_ratio))))
)

ggplot(df, aes(x,colour=Stat5_coverage)) + stat_ecdf(geom = "step") + xlab('WT TAD log(stim/naive)')
ggplot(df, aes(x,colour=Stat5_coverage)) + stat_ecdf(geom = "step") + xlab('naive TAD log(DKO/WT)')
ks.test(tad_high_ratio,tad_low_ratio)

merged_tad_tcf1_stat5 <- read.table('/media/shaoqizhu/easystore/HiC_CD8HP/TAD/merged_tad_tcf1_stat5.bed')
tad_tcf1_stat5_high <- order(merged_tad_tcf1_stat5$V7,decreasing = T)[1:200]
tad_tcf1_stat5_low <- order(merged_tad_tcf1_stat5$V7,decreasing = F)[1:200]
tad_high_ratio <- (merged_tad_pet$stim_DKO_2016[tad_tcf1_stat5_high]/merged_tad_pet$stim_WT_2016[tad_tcf1_stat5_high])
tad_low_ratio <- (merged_tad_pet$stim_DKO_2016[tad_tcf1_stat5_low]/merged_tad_pet$stim_WT_2016[tad_tcf1_stat5_low])
boxplot(tad_high_ratio,tad_low_ratio,names=c('high','low'),main='stim TAD log(dKO/WT)')

df <- data.frame(
  x = c(tad_high_ratio,tad_low_ratio),
  Tcf1_Stat5_coverage = factor(c(rep('high',length(tad_high_ratio)),rep('low',length(tad_low_ratio))))
)
ggplot(df, aes(x,colour=Tcf1_Stat5_coverage)) + stat_ecdf(geom = "step") + xlab('TAD log(dKO/WT)')
ks.test(tad_high_ratio,tad_low_ratio)



name = 'naive_WT_2019'
tad_len <- read.table(paste('/media/shaoqizhu/easystore/HiC_CD8HP/TAD/',name,'/10000_blocks.bedpe',sep = ''))
tadLen <- tad_len$V3-tad_len$V2
hist(tadLen[tadLen<1500000],main=paste(name,' TAD length',sep = ''),breaks = seq(min(tadLen),1500000,length.out=30))


naive_tad <- read.table('/media/shaoqizhu/easystore/HiC_CD8HP/TAD/naive_2016_tad.bed')
write.table(cbind(paste(naive_tad$V1,naive_tad$V2,naive_tad$V3,sep = ':'),
                  paste(naive_tad$V1,naive_tad$V2,naive_tad$V3,sep = ':')),
            '/media/shaoqizhu/easystore/HiC_CD8HP/TAD/naive_tad_straw.bed',quote = F,col.names = F,row.names = F,sep = '\t')

naive_tad_pet <- read.table('/media/shaoqizhu/easystore/HiC_CD8HP/TAD/naive_tad/naive_tad_pet.txt')
colnames(naive_tad_pet) <- c('naive_DKO_2016','naive_WT_2016')
plot(naive_tad_pet$naive_DKO_2016,naive_tad_pet$naive_WT_2016)


plot(merged_tad_pet[,1],merged_tad_pet[,3])

peak[[1]]

library(readr)
peak <- list()
peak[[1]] <-  read.table('/media/shaoqizhu/easystore/CD8-HP/CTCF/peak_diff/WT_mt_DKO_in_Tcf1_na.bed')
peak[[2]] <-  read.table('/media/shaoqizhu/easystore/CD8-HP/CTCF/peak_diff/WT_mt_DKO_no_Tcf1_na.bed')
peak[[3]] <-  read.table('/media/shaoqizhu/easystore/CD8-HP/CTCF/peak_diff/diff_kmeans_no24h/peak_diff_cluster_3.bed')
hits<- findOverlaps(GRanges(Tcf1),GRanges(merge))
peak[[4]] <-  Tcf1[-hits@from,][sample(1:nrow(Tcf1[-hits@from,]),1000),]
colnames(peak[[4]])[1:3] <- c('V1','V2','V3')

peak_pet <- list()
for (i in 1:4){
peak[[i]] <- peak[[i]][order(peak[[i]]$V1,peak[[i]]$V2),]
rownames(peak[[i]]) <- 1:nrow(peak[[i]])
peak[[i]]$V1 <- gsub('chr','',peak[[i]]$V1)
peak[[i]]$V2 <- floor(peak[[i]]$V2/10^4)*10^4
peak_pet[[i]] <- matrix(NA,nrow = nrow(peak[[i]]),ncol = 4)
}

merge <- read.csv('/media/shaoqizhu/easystore/CD8-HP/CTCF/merged_peak/merged.bed',sep='\t',header = F)
colnames(merge)[1:3] <- c('chr','start','end')
merge_bin <- merge
colnames(merge_bin)[1:3] <- c('V1','V2','V3')
merge_bin$V2 <- floor(merge$start/10^4)*10^4
merge_bin$V1 <- gsub('chr','',merge_bin$V1)

merge_bin$V2 <- merge_bin$V2 + 50000
for (i in 1:4){peak[[i]]$V2 <- peak[[i]]$V2 +50000}

mm9_size_bin <- mm9_size
mm9_size_bin$V2 <- floor(mm9_size_bin$V2/10^4)*10^4
mm9_size_bin$V1 <- gsub('chr','',mm9_size_bin$V1)

for (i in 1:4){
  for (j in unique(peak[[i]]$V1)){
    peak[[i]]$V2[peak[[i]]$V1==j][peak[[i]]$V2[peak[[i]]$V1==j]>mm9_size_bin$V2[mm9_size_bin$V1==j]] <- mm9_size_bin$V2[mm9_size_bin$V1==j]
  }
}
for (j in unique(merge_bin$V1)){
  merge_bin$V2[merge_bin$V1==j][merge_bin$V2[merge_bin$V1==j]>mm9_size_bin$V2[mm9_size_bin$V1==j]] <- mm9_size_bin$V2[mm9_size_bin$V1==j]
}

name <- c("naive_DKO_2016","naive_WT_2016","stim_DKO_2016","stim_WT_2016")

for (n in 1:length(name)){
  for (i in unique(peak[[1]]$V1)){
    hic <- read_table2(paste('/media/shaoqizhu/easystore/HiC_CD8HP/straw_KR/',name[n],'/',i,'.txt',sep=''),col_names = F)
    for (k in 1:4){
    peak_chr <- peak[[k]][peak[[k]]$V1==i,]
    for (j in 1:nrow(peak_chr)){
      #idx <- which((hic$X1 == peak_chr$V2[j] & hic$X2 %in% merge_bin$V2[merge_bin$V1==i])|
      #               (hic$X2 == peak_chr$V2[j] & hic$X1 %in% merge_bin$V2[merge_bin$V1==i]))
      idx <- which(hic$X1 == peak_chr$V2[j] | hic$X2 == peak_chr$V2[j])
      idx <- setdiff(idx,which(hic$X1 == peak_chr$V2[j] & hic$X2 == peak_chr$V2[j]))
      if(length(idx)>0){
      idx <- idx[which(abs(hic[idx,1]-hic[idx,2])<=500000)]
      if(length(idx)==0){peak_pet[[k]][as.numeric(rownames(peak_chr)[j]),n]=0}
      else{peak_pet[[k]][as.numeric(rownames(peak_chr)[j]),n] = sum(hic[idx,'X3'][hic[idx,'X3']>10],na.rm = T)}
      }else{peak_pet[[k]][as.numeric(rownames(peak_chr)[j]),n] = 0}
    }
    }
  }
}

peak_pet_ramdom <- peak_pet
peak_pet_ctcf <- peak_pet
boxplot(peak_pet_ctcf[[1]][,2],peak_pet_ramdom[[1]][,2],main='WT_mt_DKO_no_Tcf1 WT_na_2016',names=c('CTCF peaks','random peaks'))
boxplot(peak_pet_ctcf[[1]][,1],peak_pet_ctcf[[2]][,1],main='WT_na_2016',names=c('with Tcf1','without Tcf1'))
boxplot(peak_pet_ctcf[[1]][,1]-peak_pet_ctcf[[1]][,2],main='DKO_na-WT_na',names="")

peak_chr$V2[2]
idx <- which((hic$X1 == peak_chr$V2[2] & hic$X2 %in% merge_bin$V2[merge_bin$V1==i])|
               (hic$X2 == peak_chr$V2[2] & hic$X1 %in% merge_bin$V2[merge_bin$V1==i]))
as.data.frame(hic[idx,])

boxplot(log(peak_pet_ctcf[[1]][,1]/peak_pet_ctcf[[1]][,2]),main='HiC DKO-WT stim 2016 in WT_mt_DKO_no_Tcf1_na')
boxplot(peak_pet[[1]][,2],peak_pet[[4]][,2],main='HiC DKO-WT stim 2016 in WT_mt_DKO_no_Tcf1_na')

peak_hic <- list()
peak_hic$ctcf_WT_mt_DKO_na_in_Tcf1<- peak_pet
peak_hic$ctcf_WT_mt_DKO_na_no_Tcf1<- peak_pet
peak_hic$ctcf_WT_72h_high <- peak_pet

diff1 <- peak_hic$ctcf_WT_mt_DKO_na_in_Tcf1$stim_DKO_2019-peak_hic$ctcf_WT_mt_DKO_na_in_Tcf1$stim_WT_2019
diff2 <- peak_hic$ctcf_WT_mt_DKO_na_no_Tcf1$stim_DKO_2019-peak_hic$ctcf_WT_mt_DKO_na_no_Tcf1$stim_WT_2019
boxplot(diff1,diff2[diff2>-400],names=c('with Tcf1','no Tcf1'),main = 'HiC DKO-WT stim 2019 in ctcf_WT_mt_DKO_na')

boxplot(peak_hic$ctcf_WT_72h_high$stim_WT_2019-peak_hic$ctcf_WT_72h_high$naive_WT_2019
        ,main = 'HiC stim-naive WT 2016 in ctcf_WT_72h_high')

colnames(peak)[1:3] <- c('chr','start','end')

for (i in 1:7){
  hits <- findOverlaps(GRanges(peak),GRanges(diff_peak_kmeans[diff_peak_kmeans$cluster_idx==i,]))
  #print(length(unique(hits@from)))
  #print(length(unique(hits@from))/hits@nRnode)
  print(hits@nRnode)
}


