library(Seurat)
library(pheatmap)
library(stringr)
library(ggplot2)
library(openxlsx)
library(ggrepel)
library(pbapply)

Kat2a.data = Read10X('/media/shaoqizhu/easystore/Kat2a/')
cell_ko=Kat2a.data@Dimnames[[2]][grep('1',Kat2a.data@Dimnames[[2]])]
cell_wt=Kat2a.data@Dimnames[[2]][grep('2',Kat2a.data@Dimnames[[2]])]

dim(Kat2a.data[,cell_ko])
Kat2a <- CreateSeuratObject(counts = Kat2a.data, 
                              project = "10X_HdacD0")
FeatureScatter(object = Kat2a, 
               feature1 = 'nCount_RNA', 
               feature2 = 'nFeature_RNA')

Kat2a[["percent.mt"]] <- PercentageFeatureSet(Kat2a, pattern = "^mt-")

VlnPlot(Kat2a, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

Kat2a <- subset(Kat2a, subset = percent.mt < 5)

Kat2a <- NormalizeData(object = Kat2a, 
                         normalization.method = "LogNormalize",
                         scale.factor = 10000)
Kat2a <- FindVariableFeatures(object = Kat2a, 
                                mean.function = ExpMean, 
                                dispersion.function = LogVMR, 
                                x.low.cutoff = 0, 
                                x.high.cutoff = 4, 
                                y.cutoff = 0.5)
VariableFeaturePlot(object = Kat2a)
all.genes<-rownames(Kat2a)
Kat2a <- ScaleData(object = Kat2a, features = all.genes)
Kat2a <- RunPCA(object = Kat2a, features = VariableFeatures(object = Kat2a))
print(x = Kat2a[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object = Kat2a, dims = 1:2, reduction = "pca")
Kat2a <- FindNeighbors(object = Kat2a, dims = 1:10)
Kat2a <- RunTSNE(object = Kat2a, dims.use = 1:10)
Kat2a <- FindClusters(object = Kat2a,resolution = 0.25)
DimPlot(Kat2a,reduction="tsne",label=T,pt.size=0.5,cells.highlight = group2) 
FeaturePlot(Kat2a,'Rps7')

hist(Kat2a@assays$RNA@counts['Kat2a',colnames(Kat2a@assays$RNA@counts)%in%cell_ko],main='group 1')
hist(Kat2a@assays$RNA@counts['Kat2a',colnames(Kat2a@assays$RNA@counts)%in%cell_wt],main='group 2')
length(group1)

mtx_ket2a <- Kat2a@assays$RNA@counts[,colnames(Kat2a@assays$RNA@counts)%in%cell_wt]
mtx_ket2a <- t(t(mtx_ket2a)/apply(mtx_ket2a, 2, sum)*mean(apply(mtx_ket2a, 2, sum)))

cov <- matrix(data=0,nrow = nrow(mtx_ket2a),ncol=2)
cov_2states <- pbsapply(
  X = c('Kat2a'),FUN = function(x) {
    cell_top <- names(sort(mtx_ket2a[x,],decreasing = T)[1:100])
    cell_bot <- names(sort(mtx_ket2a[x,],decreasing = F)[1:100])
    mean(mtx_ket2a[x,cell_top])
    cov[,1] <- apply(mtx_ket2a[,cell_bot],1,function(x) sd(x)/mean(x))
    cov[,2] <- apply(mtx_ket2a[,cell_top],1,function(x) sd(x)/mean(x))
    c(mean(mtx_ket2a[x,cell_bot]),mean(mtx_ket2a[x,cell_top]),
      length(which(cov[,1] > cov[,2])),length(which(cov[,1] < cov[,2])))}
)
t(cov_2states)

counts <- as.matrix(Kat2a@assays$RNA@counts)
counts <- t(t(counts)/apply(counts, 2, sum)*mean(apply(counts, 2, sum)))
cov_wt <- apply(counts[,cell_wt]+0.1,1,function(x) sd(x)/mean(x))
cov_ko <- apply(counts[,cell_ko]+0.1,1,function(x) sd(x)/mean(x))
ratio <- log((cov_ko+0.1)/(cov_wt+0.1))
hist(ratio[ratio>-1&ratio<1&ratio!=0],main='hist of log(Kat2a_ko/Kat2a_wt) (sd/mean)',xlab='log(ratio)')

x = c('Kat2a')
  cell_top <- names(sort(mtx_ket2a[x,],decreasing = T)[1:100])
  cell_bot <- names(sort(mtx_ket2a[x,],decreasing = F)[1:100])
  mean(mtx_ket2a[x,cell_top])
  cov[,1] <- apply(mtx_ket2a[,cell_bot],1,function(x) sd(x)/mean(x))
  cov[,2] <- apply(mtx_ket2a[,cell_top],1,function(x) sd(x)/mean(x))
  c(mean(mtx_ket2a[x,cell_bot]),mean(mtx_ket2a[x,cell_top]),
    length(which(cov[,1] > cov[,2])),length(which(cov[,1] < cov[,2])))

A <- length(which(cov[,1] > cov[,2]));B <- length(which(cov[,1] < cov[,2]));
a <- length(which(ratio>0));b <- length(which(ratio<0))
Aa <- length(intersect(rownames(mtx_ket2a)[which(cov[,1] > cov[,2])],names(ratio)[ratio>0]))
Ab <- length(intersect(rownames(mtx_ket2a)[which(cov[,1] > cov[,2])],names(ratio)[ratio<0]))
Ba <- length(intersect(rownames(mtx_ket2a)[which(cov[,1] < cov[,2])],names(ratio)[ratio>0]))
Bb <- length(intersect(rownames(mtx_ket2a)[which(cov[,1] < cov[,2])],names(ratio)[ratio<0]))

cbind(c(Aa,Ab,A),c(Ba,Bb,B),c(a,b,NA))

fisher.test(cbind(c(5859,545),c(1417,191)))





