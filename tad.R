library(readr)
library(pheatmap)
merged_tad <- read.table('/media/shaoqizhu/easystore/HiC_CD8HP/TAD_naive/WT_CD8_naive_pooled_TAD.bedpe')
merged_tad <- merged_tad[,1:3]
name <- c("naive_DKO_2016","naive_DKO_2019","naive_WT_2016","naive_WT_2019")
intra_tad <- matrix(NA,nrow = nrow(merged_tad),ncol = 4)
total_tad <- matrix(NA,nrow = nrow(merged_tad),ncol = 4)

for (n in 1:length(name)){
  for (i in unique(merged_tad$V1)){
    hic <- read_table2(paste('/media/shaoqizhu/easystore/HiC_CD8HP/hic_straw/',name[n],'/',i,'.txt',sep=''),col_names = F)
    merged_tad_chr <- merged_tad[merged_tad$V1==i,]
    for (j in 1:nrow(merged_tad_chr)){
      intra_tad_idx <- which((hic$X2 >= merged_tad_chr$V2[j] & hic$X2 <= merged_tad_chr$V3[j]) & 
                         (hic$X3 >= merged_tad_chr$V2[j] & hic$X3 <= merged_tad_chr$V3[j]))
      intra_tad_idx <- setdiff(intra_tad_idx,which(hic$X2==hic$X3))
      
      total_tad_idx <- which((hic$X2 >= merged_tad_chr$V2[j] & hic$X2 <= merged_tad_chr$V3[j]) | 
                         (hic$X3 >= merged_tad_chr$V2[j] & hic$X3 <= merged_tad_chr$V3[j]))
      total_tad_idx <- setdiff(total_tad_idx,which(hic$X2==hic$X3))
      
      intra_tad[rownames(merged_tad)==rownames(merged_tad_chr)[j],n] = sum(hic[intra_tad_idx,'X4'])
      total_tad[rownames(merged_tad)==rownames(merged_tad_chr)[j],n] = sum(hic[total_tad_idx,'X4'])
    }
  }
}
tad_total <- merged_tad[,4:7]

tad_idx <- which((hic$X2 >= merged_tad_chr$V2[1] & hic$X2 <= merged_tad_chr$V3[1]) | 
        (hic$X3 >= merged_tad_chr$V2[1] & hic$X3 <= merged_tad_chr$V3[1]))
hic[tad_idx,'X4']
merged_tad_chr[1,]
#=================================== random tad
merged_tad <- read.table('/media/shaoqizhu/easystore/HiC_CD8HP/TAD/merged_tad.bed')
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

random_tad <- as.data.frame(random_tad)
name <- c("naive_DKO_2016","naive_DKO_2019","naive_WT_2016","naive_WT_2019","stim_DKO_2016","stim_DKO_2019","stim_WT_2016","stim_WT_2019")
for (n in 1:length(name)){
  for (i in c(1:19,'X')){
    hic <- read_table2(paste('/media/shaoqizhu/easystore/HiC_CD8HP/hic_straw/',name[n],'/',i,'.txt',sep=''),col_names = F)
    random_tad_chr <- random_tad[random_chr==i,]
    for (j in 1:nrow(random_tad_chr)){
      tad_idx <- which((hic$X2 >= random_tad_chr$V1[j] & hic$X2 <= random_tad_chr$V2[j]) | 
                         (hic$X3 >= random_tad_chr$V1[j] & hic$X3 <= random_tad_chr$V2[j]))
      if (length(tad_idx)==0){random_tad[rownames(random_tad)==rownames(random_tad_chr)[j],n+2]=0}
      else {random_tad[rownames(random_tad)==rownames(random_tad_chr)[j],n+2] = sum(hic[tad_idx,'X4'])}
    }
  }
}

merged_tad_pet <- read.table('/media/shaoqizhu/easystore/HiC_CD8HP/TAD_naive/merged_tad_naive.bed')
colnames(merged_tad_pet)[1:4] <- c("naive_DKO_2016","naive_DKO_2019","naive_WT_2016","naive_WT_2019")

merged_tad_pet <- merged_tad_pet/tad_total
merged_tad_pet <- as.data.frame(t(t(merged_tad_pet)/apply(merged_tad_pet, 2, sum)*10^4))
merged_tad_pet_Z <- (merged_tad_pet-apply(X = merged_tad_pet, MARGIN = 1, mean))/apply(X = merged_tad_pet, MARGIN = 1, sd)

set.seed(2);TAD_kmeans <- kmeans(merged_tad_pet_Z,centers=5,nstart=25)
TAD_kmeans$reorder <- TAD_kmeans$cluster; reorders <- c(1,2,3,4,5)# 2019 reorders <- c(5,1,4,2,3);  2016 reorders <- c(4,2,3,5,1)
for (i in 1:5){TAD_kmeans$reorder[TAD_kmeans$cluster==i] <- reorders[i]}
tad_kmeans <- as.data.frame(sort(TAD_kmeans$reorder)); names(tad_kmeans) <- 'cluster'; 
tad_kmeans$cluster <- factor(tad_kmeans$cluster)
rownames(tad_kmeans) <- order(TAD_kmeans$reorder)
pheatmap(merged_tad_pet_Z[order(TAD_kmeans$reorder),], cluster_rows = F,cluster_cols = F,show_rownames = F,
         angle_col='315',annotation_row = tad_kmeans)
