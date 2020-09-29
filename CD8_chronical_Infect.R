library(Seurat)
library(pheatmap)
library(stringr)
library(ggplot2)
library(openxlsx)
library(ggrepel)
library(pbapply)

CD8.data = Read10X('/media/shaoqizhu/easystore/CD8_chronical_infection/GSE122675/')

CD8 <- CreateSeuratObject(counts = CD8.data, 
                            project = "10X_HdacD0")
FeatureScatter(object = CD8, 
               feature1 = 'nCount_RNA', 
               feature2 = 'nFeature_RNA')

CD8[["percent.mt"]] <- PercentageFeatureSet(CD8, pattern = "^mt-")

VlnPlot(CD8, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#CD8 <- subset(CD8, subset = percent.mt < 5)

CD8 <- NormalizeData(object = CD8, 
                       normalization.method = "LogNormalize",
                       scale.factor = 10000)
CD8 <- FindVariableFeatures(object = CD8, 
                              mean.function = ExpMean, 
                              dispersion.function = LogVMR, 
                              x.low.cutoff = 0, 
                              x.high.cutoff = 4, 
                              y.cutoff = 0.5)
VariableFeaturePlot(object = CD8)
all.genes<-rownames(CD8)
CD8 <- ScaleData(object = CD8, features = all.genes)
CD8 <- RunPCA(object = CD8, features = VariableFeatures(object = CD8))
print(x = CD8[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object = CD8, dims = 1:2, reduction = "pca")
CD8 <- FindNeighbors(object = CD8, dims = 1:10)
CD8 <- RunTSNE(object = CD8, dims.use = 1:10)
CD8 <- FindClusters(object = CD8,resolution = 0.1)
DimPlot(CD8,reduction="tsne",label=T,pt.size=0.5,cells.highlight = names(which(CD8@assays$RNA@counts['Tcf7',]>1))) 
FeaturePlot(CD8,features = 'Tcf7')
FeaturePlot(CD8,features = 'Entpd1')
table(CD8@assays$RNA@counts['Entpd1',which(CD8@assays$RNA@counts['Tcf7',]>=1)])
hist(CD8@assays$RNA@counts['Entpd1',which(CD8@assays$RNA@counts['Tcf7',]>=1)],xlab='Entpd1',main='Entpd1 in Tcf7 cluster')
CD8@assays$RNA@counts[]
names(which(CD8@assays$RNA@counts['Tcf7',]>1))
names(CD8@active.ident)[CD8@active.ident==2]]