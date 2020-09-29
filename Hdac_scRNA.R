library(Seurat)
library(pheatmap)
library(dplyr)
library(stringr)
library(devtools)
library(ggplot2)
library(openxlsx)
library(ggrepel)
library(cowplot)
library(pbapply)
fileName_D0 = "/media/shaoqizhu/easystore/Hdac-scRNAseq/Ctrl-D0_20180801002_filtered_gene.h5"
fileName_D1 = "/media/shaoqizhu/easystore/Hdac-scRNAseq/Ctrl-D1_20180801000_filtered_gene.h5"
fileName_D0_Dko = "/media/shaoqizhu/easystore/Hdac-scRNAseq/Dko-D0_20180801000_filtered_gene.h5"
fileName_D1_Dko = "/media/shaoqizhu/easystore/Hdac-scRNAseq/Dko-D1_20180801002_filtered_gene.h5"

HdacD0.data <- Read10X_h5(fileName_D0)
HdacD0_Dko.data <- Read10X_h5(fileName_D0_Dko)
HdacD1.data <- Read10X_h5(fileName_D1)
HdacD1_Dko.data <- Read10X_h5(fileName_D1_Dko)

HdacD0.data@Dimnames[[2]] <- paste(str_split_fixed(HdacD0.data@Dimnames[[2]],'-',2)[,1],'1',sep='-')
HdacD0_Dko.data@Dimnames[[2]] <- paste(str_split_fixed(HdacD0_Dko.data@Dimnames[[2]],'-',2)[,1],'2',sep='-')
HdacD1.data@Dimnames[[2]] <- paste(str_split_fixed(HdacD1.data@Dimnames[[2]],'-',2)[,1],'3',sep='-')
HdacD1_Dko.data@Dimnames[[2]] <- paste(str_split_fixed(HdacD1_Dko.data@Dimnames[[2]],'-',2)[,1],'4',sep='-')

HdacD01.data <- as(cbind(as.matrix(HdacD0.data),as.matrix(HdacD0_Dko.data),as.matrix(HdacD1.data),as.matrix(HdacD1_Dko.data)), "dgCMatrix")

HdacD01 <- CreateSeuratObject(counts = HdacD01.data, 
                              min.cells = 100, min.features = 500,
                              project = "10X_HdacD0")
FeatureScatter(object = HdacD01, 
               feature1 = 'nCount_RNA', 
               feature2 = 'nFeature_RNA')

HdacD01[["percent.mt"]] <- PercentageFeatureSet(HdacD01, pattern = "^mt-")

VlnPlot(HdacD01, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

HdacD01 <- subset(HdacD01, subset = percent.mt < 5)

HdacD01 <- NormalizeData(object = HdacD01, 
                         normalization.method = "LogNormalize",
                         scale.factor = 10000)
HdacD01 <- FindVariableFeatures(object = HdacD01, 
                                mean.function = ExpMean, 
                                dispersion.function = LogVMR, 
                                x.low.cutoff = 0, 
                                x.high.cutoff = 4, 
                                y.cutoff = 0.5)
VariableFeaturePlot(object = HdacD01)
all.genes<-rownames(HdacD01)
HdacD01 <- ScaleData(object = HdacD01, features = all.genes)
HdacD01 <- RunPCA(object = HdacD01, features = VariableFeatures(object = HdacD01))
print(x = HdacD01[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object = HdacD01, dims = 1:2, reduction = "pca")
HdacD01 <- FindNeighbors(object = HdacD01, dims = 1:10)
HdacD01 <- RunTSNE(object = HdacD01, dims.use = 1:10)
HdacD01 <- FindClusters(object = HdacD01,resolution = 0.25)
DimPlot(HdacD01,reduction="tsne",label=F,pt.size=0.5) 

condition <- c('WT_D0','dKO_D0','WT_D1','dKO_D1')
tsne <- as.data.frame(HdacD01@reductions$tsne@cell.embeddings); tsne$Condition <-0
for (i in 1:4){tsne$Condition[grep(i,rownames(tsne))] <- condition[i]}; tsne$Condition <- factor(tsne$Condition)
ggplot(tsne,aes(x=tSNE_1,y=tSNE_2)) + geom_point(size=0.5,aes(color=Condition)) + labs(x = "tSNE_1",y="tSNE_2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

count_WT_D0 <- as.matrix(HdacD01@assays$RNA@counts[,tsne$Condition%in%c('WT_D1','dKO_D1')])
count_WT_D0 <- rbind(GeneID=c(rep('WT_D1',671),rep('dKO_D1',938)),count_WT_D0)
write.table(count_WT_D0,'/media/shaoqizhu/easystore/Burst/Hdac_D1.txt',sep = '\t',quote = F,col.names = F)

dim(HdacD01@assays$RNA@counts[,tsne$Condition=='WT_D0'])

HdacD01.marker <- FindAllMarkers(object = HdacD01,only.pos = T)
TF <- read.csv('~/Documents/Mus_musculus_TF.txt',sep = '\t')
TF_marker <- intersect(HdacD01.marker$gene, TF$Symbol)
write.csv(HdacD01.marker,'/media/shaoqizhu/easystore/Hdac-scRNAseq/marker.csv')

pheatmap(HdacD01@assays$RNA@data[TF_marker,],cluster_rows=T,cluster_cols = F,show_rownames = T,show_colnames = F,
         main = 'TF in all cells')

FeaturePlot(HdacD01,transporter[1:12],reduction="tsne",cols =c('lightgrey','red'),ncol = 4)

top10 <- HdacD01.marker %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(HdacD01, features = top10$gene) + NoLegend()


HdacD01.marker[HdacD01.marker$gene=='Slc3a2',]
transporter <- HdacD01.marker$gene[grep('Slc',HdacD01.marker$gene)]

gene_filtered <- setdiff(HdacD01.data@Dimnames[[1]],rownames(HdacD01@assays$RNA@counts))
gene_filtered_sum <- apply(HdacD01.data[gene_filtered,],1,function(x) sum(x>0))
hist(gene_filtered_sum)
enrichr <- read.csv('/media/shaoqizhu/easystore/Hdac-scRNAseq/ChEA_2016_table (2).txt',sep='\t')
grep('NFKB',enrichr$Term)

expr_mtx <- HdacD01@assays$RNA@counts[,tsne$Condition=='WT_D0']
expr_mtx <- t(t(expr_mtx)/apply(expr_mtx, 2, sum)*10^4)
cov <- matrix(data=0,nrow = nrow(expr_mtx),ncol=2)
cov_2states <- pbsapply(
X = 'Hdac1',FUN = function(x) {
cell_top <- names(sort(expr_mtx[x,],decreasing = T)[1:100])
cell_bot <- names(sort(expr_mtx[x,],decreasing = F)[1:100])
mean(expr_mtx[x,cell_top])
cov[,1] <- apply(expr_mtx[,cell_bot],1,function(x) sd(x)/mean(x))
cov[,2] <- apply(expr_mtx[,cell_top],1,function(x) sd(x)/mean(x))
c(mean(expr_mtx[x,cell_bot]),mean(expr_mtx[x,cell_top]),
  length(which(cov[,1] > cov[,2])),length(which(cov[,1] < cov[,2])))}
)
t(cov_2states)

cov_2state <- t(cov_2states)
rownames(cov_2state) <- rownames(expr_mtx)
cov_2state[order(cov_2state[,3],decreasing = T)[1:100],]
cell_top <- names(sort(expr_mtx['Hdac1',],decreasing = T)[1:100])
cell_bot <- names(sort(expr_mtx['Hdac1',],decreasing = F)[1:100])
write.csv(cov_2state[order(cov_2state[,3],decreasing = F)[1:100],],'~/Documents/gene_cov.csv')

counts <- as.matrix(HdacD01@assays$RNA@counts)
counts <- t(t(counts)/apply(counts, 2, sum)*10^4)
cov_WT_D0 <- apply(counts[,tsne$Condition=='WT_D0'],1,function(x) sd(x)/mean(x))
cov_WT_D1 <- apply(counts[,tsne$Condition=='WT_D1'],1,function(x) sd(x)/mean(x))
cov_dKO_D0 <- apply(counts[,tsne$Condition=='dKO_D0'],1,function(x) sd(x)/mean(x))
cov_dKO_D1 <- apply(counts[,tsne$Condition=='dKO_D1'],1,function(x) sd(x)/mean(x))
mean_WT_D0 <- apply(counts[,tsne$Condition=='WT_D0'],1,function(x) mean(x))
mean_WT_D1 <- apply(counts[,tsne$Condition=='WT_D1'],1,function(x) mean(x))
mean_dKO_D0 <- apply(counts[,tsne$Condition=='dKO_D0'],1,function(x) mean(x))
mean_dKO_D1 <- apply(counts[,tsne$Condition=='dKO_D1'],1,function(x) mean(x))

hist(mean_WT_D0[mean_WT_D0<1])
ratio <- log(cov_WT_D1/cov_WT_D0)
hist(ratio[ratio>-1&ratio<1&mean_WT_D1>0.5&mean_WT_D0],breaks = seq(-1,1,0.1),
     main='hist of log(WT_D1/WT_D0) (sd/mean)',xlab='log(ratio)')
intersect(names(ratio),rownames(expr_mtx))
ratio0 <- log(cov_dKO_D0/cov_WT_D0);ratio0 <- ratio0[order(ratio0,decreasing = T)]
ratio1 <- log(cov_dKO_D1/cov_WT_D1);ratio1 <- ratio1[order(ratio1,decreasing = T)]
gene_cov <- intersect(names(ratio0[1:7395]),names(ratio1[1:2761]))

intersect(gene_cov,HdacD01.marker$gene[HdacD01.marker$cluster==0])
length(HdacD01.marker$gene[HdacD01.marker$cluster==0])

cell_bot <- sample(intersect(names(sort(expr_mtx['Hdac1',],decreasing = F)[1:300]),
                      names(sort(expr_mtx["Hdac2",],decreasing = F)[1:300])),100)
cell_top <- sample(intersect(names(sort(expr_mtx['Hdac1',],decreasing = T)[1:300]),
                      names(sort(expr_mtx["Hdac2",],decreasing = T)[1:300])),100)
mean(expr_mtx['Hdac2',cell_top])

fisher.test(cbind(c(6387,657),c(929,176)))
chisq.test(cbind(c(6387,657),c(929,176)))

ratio <- ratio[order(names(ratio))]
expr_mtx <- expr_mtx[order(rownames(expr_mtx)),]
cov[,1] <- apply(expr_mtx[,cell_bot],1,function(x) sd(x)/mean(x))
cov[,2] <- apply(expr_mtx[,cell_top],1,function(x) sd(x)/mean(x))
c(mean(expr_mtx['Hdac1',cell_bot]),mean(expr_mtx['Hdac1',cell_top]),
  length(which(cov[,1] > cov[,2])),length(which(cov[,1] < cov[,2])))

A <- length(which(cov[,1] > cov[,2]));B <- length(which(cov[,1] < cov[,2]));
a <- length(which(ratio>0));b <- length(which(ratio<0))
Aa <- length(intersect(rownames(expr_mtx)[which(cov[,1] > cov[,2])],names(ratio)[ratio>0]))
Ab <- length(intersect(rownames(expr_mtx)[which(cov[,1] > cov[,2])],names(ratio)[ratio<0]))
Ba <- length(intersect(rownames(expr_mtx)[which(cov[,1] < cov[,2])],names(ratio)[ratio>0]))
Bb <- length(intersect(rownames(expr_mtx)[which(cov[,1] < cov[,2])],names(ratio)[ratio<0]))

cbind(c(Aa,Ab,A),c(Ba,Bb,B),c(a,b,NA))

intersect(names(sort(ratio,decreasing = T)[1:1000]),rownames(expr_mtx)[order(cov[,1],decreasing = T)[1:1000]])

sd(expr_mtx[order(cov[,1],decreasing = T)[1],cell_bot])/mean(expr_mtx[order(cov[,1],decreasing = T)[1],cell_bot])
sort(cov[,1],decreasing = T)[1]

gene.names <- intersect(rownames(expr_mtx)[which(cov[,1] > cov[,2])],names(ratio)[ratio>0])

Hdac_cov <- (cov[,1]/cov[,2])[rownames(expr_mtx)%in%gene.names]
ratio_cov <- ratio[names(ratio)%in%gene.names]
which(Hdac_cov>10)
plot(log(Hdac_cov[-2297]),ratio_cov[-2297],main=paste('cor',cov(Hdac_cov[-2297],ratio_cov[-2297])),
     xlab='Hdac bot/top cov (log)',ylab='dKO_D0/WT_D0 cov (log)')
cov(Hdac_cov[-2297],ratio_cov[-2297])



{File = '3'
  GO <- read.csv(paste('/media/shaoqizhu/easystore/Hdac-scRNAseq/DAVID/',File,'.txt',sep=''),header = T, sep = "\t")
  GO$Term <- str_split_fixed(GO$Term,'~',2)[,2]
  GO$Term <- factor(GO$Term, levels=rev(unique(GO$Term)))
  GO <-GO[1:40,]
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

gene_info <- read.csv('/media/shaoqizhu/easystore/Hdac-scRNAseq/DAVID/gene_info.txt',sep='\t',header = F)

burst_D1 <- read.table('/media/shaoqizhu/easystore/Burst/Hdac/WT_D1_dKO_D1.out',sep = '\t',header=T,row.names = 1)
burst_D0 <- read.table('/media/shaoqizhu/easystore/Burst/Hdac/WT_D0_dKO_D0.out',sep = '\t',header=T,row.names = 1)
intersect(rownames(burst_D1)[burst_D1$Rf>(1)],)
plot(log(burst_D1$f1),log(burst_D1$f2))
hist(log(burst_D0$cv2/burst_D0$cv1))
hist(log(burst_D1$cv2/burst_D1$cv1))
marker3 <- HdacD01.marker$gene[HdacD01.marker$cluster==3]
gene = rownames(burst_D1)[order(burst_D1$s1,decreasing = F)]
expr <- HdacD01@assays$RNA@counts[gene[1],tsne$Condition=='dKO_D1']
hist(expr,breaks = seq(0,max(expr),1),main=gene[9])
RidgePlot(HdacD01,features = gene[1:9], idents=3, slot = 'counts',same.y.lims=T)
sort(burst_D1$s1,decreasing = F)[1:9]
hist(burst_D1$Rs[rownames(burst_D1)%in%marker3])
hist(burst_D1$Rf[rownames(burst_D1)%in%marker3])
