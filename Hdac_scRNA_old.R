library(rhdf5)
library(Seurat)
library(pheatmap)
library(dplyr)
library(stringr)
library(devtools)
library(ggplot2)
library(openxlsx)
library(ggrepel)
library(corrplot)
library(seriation)
library(cowplot)
library(pbapply)
fileName_D0 = "/media/shaoqizhu/easystore/hdac_scRNA_seq/Ctrl-D0_20180801002_filtered_gene.h5"
fileName_D1 = "/media/shaoqizhu/easystore/hdac_scRNA_seq/Ctrl-D1_20180801000_filtered_gene.h5"
fileName_D0_Dko = "/media/shaoqizhu/easystore/hdac_scRNA_seq/Dko-D0_20180801000_filtered_gene.h5"
fileName_D1_Dko = "/media/shaoqizhu/easystore/hdac_scRNA_seq/Dko-D1_20180801002_filtered_gene.h5"

data_D0_Dko <- h5read(fileName_D0_Dko,'mm10')
data1 <- h5read(fileName1,'matrix')

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
HdacD01 <- FindClusters(object = HdacD01,resolution = 0.9)
DimPlot(HdacD01,reduction="tsne",label=T) 
HdacD01 <- RunUMAP(object = HdacD01, dims = 1:10)

DimPlot(HdacD01,reduction="tsne",cols.highlight=c("grey","red","blue","yellow"),
        cells.highlight = list(HdacD1_Dko.data@Dimnames[[2]],
                               HdacD0_Dko.data@Dimnames[[2]],
                               HdacD1.data@Dimnames[[2]],
                               HdacD0.data@Dimnames[[2]]))

HdacD01.DEG <- FindAllMarkers(object = HdacD01,only.pos = T)

top10 <- HdacD01.DEG %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(HdacD01, cells=cluster5_wtD0,#names(HdacD01@active.ident[HdacD01@active.ident%in%c(0,2,5)]),
          features = HdacD01.DEG$gene[HdacD01.DEG$cluster==0])

col_group <-data.frame(group=c(rep('wtD0',length(cluster5_wtD0)),rep('dkoD0',length(cluster5_DkoD0))))
rownames(col_group) <- c(cluster5_wtD0,cluster5_DkoD0)
mat <- as.data.frame(HdacD01@assays$RNA@data[HdacD01.DEG$gene[HdacD01.DEG$cluster==0],c(cluster5_wtD0,cluster5_DkoD0)])
mat_mean = apply(X = mat, MARGIN = 1, mean)
mat_zscore <- mat/mat_mean-1
mat_zscore[mat_zscore>10] <- 0

pheatmap(mat_zscore,cluster_rows=F,cluster_cols=F,show_rownames=T,show_colnames=F,annotation_col = col_group)

table(HdacD01.DEG$cluster)

cluster5_wtD0 <- intersect(names(HdacD01@active.ident[HdacD01@active.ident==5]),HdacD0.data@Dimnames[[2]])
cluster5_DkoD0 <- intersect(names(HdacD01@active.ident[HdacD01@active.ident==5]),HdacD0_Dko.data@Dimnames[[2]])
sum(HdacD01@assays$RNA@data['Hdac1',cluster5_wtD0]>0.5)/length(HdacD01@assays$RNA@data['Hdac1',cluster5_wtD0])
sum(HdacD01@assays$RNA@data['Hdac1',cluster5_DkoD0]>0.5)/length(HdacD01@assays$RNA@data['Hdac1',cluster5_DkoD0])


FeaturePlot(HdacD01,features='Hdac1',reduction="tsne")
FeaturePlot(HdacD01,features='Hdac2',reduction="tsne")

expMtx <- as.data.frame(HdacD01@assays$RNA@data)
marker = data.frame(avg_logFC=numeric(),p_val=numeric(),
                    p_val_adj=numeric(),cluster = numeric(),
                    gene = character())
for(i in 1:9){
  thres = 0.4
  idx = as.numeric(HdacD01@active.ident)
  cell.1 = names(expMtx)[idx==i]
  cell.2 = names(expMtx)[idx!=i]
  data.1 = apply(X = expMtx[, cell.1, drop = F], MARGIN = 1, 
                 FUN = function(x) log2(x = mean(x) + 1))
  data.2 = apply(X = expMtx[, cell.2, drop = F], MARGIN = 1, 
                 FUN = function(x) log2(x = mean(x) + 1))
  avg_logFC = data.1-data.2
  genes = names(avg_logFC)[avg_logFC>thres]
  group.info <- data.frame(row.names = c(cell.1, cell.2))
  group.info[cell.1, "group"] <- "Group1"
  group.info[cell.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(x = group.info[, "group"])
  p_val <- pbsapply(
    X = genes,FUN = function(x) {
      return(wilcox.test(t(expMtx[x,rownames(group.info)]) ~ group.info[, "group"])$p.value)})
  p_val_adj <- p.adjust(p = p_val,method = "bonferroni",n = length(genes))
  
  DEG = as.data.frame(cbind(avg_logFC = signif(avg_logFC[genes],7),
                            p_val = signif(unlist(p_val),7),
                            p_val_adj=signif(p_val_adj,7),
                            cluster=rep(i-1,length(genes)),
                            gene=genes))
  marker = rbind(marker,DEG)
}
table(marker$cluster)
write.csv(gene_module,'~/Downloads/module_genes_log.csv')
top10 <- marker %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(HdacD01, cells=names(expMtx)[as.numeric(HdacD01@active.ident)%in%c(0,2,5)],
          features = as.character(marker$gene)[marker$cluster%in%c(0,2,5)]) + NoLegend()



ExpMtxAll <- as.data.frame(t(HdacD01@assays$RNA@data[unique(HdacD01.DEG$gene),]))
ExpMtxAll$ident <- as.numeric(HdacD01@active.ident)
ExpMean <- aggregate(ExpMtxAll,list(ExpMtxAll$ident),mean)[,c(-1,-length(ExpMtxAll)-1)]

#=================Fisher test
p_val_fisher <- list()
for (i in 0:13){
fisherTest <- function(cell){
  geneRef <- rownames(gene_module)[gene_module$module==i]
  #as.character(marker$gene[marker$cluster==1])#
  geneNonzero = rownames(expMtx)[expMtx[,cell]>3]
  interNum <- length(intersect(geneNonzero,geneRef))
  test <- fisher.test(matrix(c(interNum, length(geneRef)-interNum,length(geneNonzero)-interNum, 
                       nrow(expMtx)-length(geneRef)-length(geneNonzero)+interNum),
                     nrow=2),alternative = 'greater')
  test$p.value
}
p_val_fisher[[i+1]] <- pbsapply(X = 1:length(expMtx),FUN = fisherTest)
}
p_val_fisher <- as.data.frame(p_val_fisher)
names(p_val_fisher) = 0:13
tsne <- DimPlot(HdacD01,reduction="tsne")$data
tsne <- cbind(tsne,p_val_fisher)

library(cowplot)
library(gridExtra)
g <- list()
for (i in 1:14){
  Data = tsne[,c(1,2,i+3)]; #names(Data)[3]<-'pvalue'
  ggplot(data = Data, aes(x = Data$tSNE_1,y = Data$tSNE_2),) +
  geom_point(data = Data,aes(color = -1*log10(Data[,3]))) +
  scale_color_gradient(low="white",high="red") +
  labs(color = expression(-log[10](p.value)),x = "tSNE_1",y="tSNE_2",
       title=paste("module", i-1 ,"genes")) +
  theme(plot.title=element_text(size=12,vjust=0.5,hjust=0.5))
  ggsave(paste('~/Document/fisher test log/module ',i-1,'.png',sep=''))
}
grid.arrange(grobs=g,ncol = 3)

Intersect <- as.data.frame(matrix(data=NA,nrow = 9,ncol = 14))
for (i in 0:8){
  for (j in 0:13){
    Intersect[i+1,j+1] <- length(intersect(as.character(HdacD01.DEG$gene[HdacD01.DEG$cluster==i]),
              rownames(gene_module)[gene_module$module==j]))
    }
}
colnames(Intersect) <- 0:13
rownames(Intersect) <- 0:8
pheatmap(Intersect,cluster_rows=F,cluster_cols = F,display_numbers=T,
         number_format='%.0f',fontsize=15)

idx = as.numeric(HdacD01@active.ident)
DimPlot(HdacD01,reduction="tsne",cells.highlight = names(expMtx)[idx==1])

#=================GO enrichment plot

for (i in c(1:2,4:9,12)){
  GO <- read.csv(paste('~/Downloads/GO1/enrichment.KEGG (',i,').tsv',sep=''),header = T, sep = "\t")
  GO$term.description <- factor(GO$term.description, levels=rev(unique(GO$term.description)))
  GO$ratio <- GO$observed.gene.count/GO$background.gene.count
  GO <-GO[1:30,]
  ggplot(data = GO, aes(x = GO$observed.gene.count,y = GO$term.description),) +
    geom_point(aes(size = -1*log10(GO$false.discovery.rate),
                   color = GO$ratio)) +
    scale_color_gradient(low="green",high="red") +
    labs(color = "Gene Ratio",size = expression(-log[10](FDR)),
         x = "Gene Count",y="GO term",
         title=paste('Module',i,'KEGG',sep=' ')) +
    theme(plot.title=element_text(size=12,vjust=0.5,hjust=0.5),
          legend.title=element_text(size=10),
          legend.text=element_text(size=10),
          legend.position="right",
          axis.text=element_text(size=10),
          axis.title=element_text(size=10),)
  ggsave(paste('~/Downloads/GO/KEGG module ',i,'.png',sep=''))
}

#Comparison of gene expression for different WGCNA modules between Ctrl1 and Dko1
ExpMtx <- as.data.frame(HdacD01@assays$RNA@data[rownames(gene_module),
                                                  c(HdacD1.data@Dimnames[[2]],HdacD1_Dko.data@Dimnames[[2]])])
gene_module$module <- as.factor(gene_module$module)
cell_sample <- data.frame(sample=rep(c('Ctrl1','Dko1'),c(length(HdacD1.data@Dimnames[[2]]),length(HdacD1_Dko.data@Dimnames[[2]]))))
cell_sample$sample <- as.factor(cell_sample$sample)
rownames(cell_sample) <- colnames(ExpMtx)
pheatmap(ExpMtx_Z,cluster_rows=F,cluster_cols=F,show_rownames=F,show_colnames=F,
         annotation_row = gene_module,annotation_col = cell_sample)

max(ExpMtx)

cell.1 = HdacD1.data@Dimnames[[2]]; cell.2 = HdacD1_Dko.data@Dimnames[[2]];
data.1 = apply(X = ExpMtx[, cell.1, drop = F], MARGIN = 1, 
               FUN = function(x) x = mean(x = expm1(x = x)))
data.2 = apply(X = ExpMtx[, cell.2, drop = F], MARGIN = 1, 
               FUN = function(x) x = mean(x = expm1(x = x)))
avg_logFC = data.1/data.2

Mean = apply(X = ExpMtx, MARGIN = 1, mean)
ExpMtx_Z <- ExpMtx/Mean-1
ExpMtx_Z[ExpMtx_Z>8] <- 0

for (i in 0:7){print(mean(avg_logFC[gene_module$module==i]))}

g <- list(); for (i in 0:7) {
avg <- avg_logFC[gene_module$module==i]
h <- hist(avg[avg<4.6],breaks = seq(0.8,4.6,0.2),plot=F); 
Data <- as.data.frame(cbind(h$mids,h$counts))
g[[i+1]] <- ggplot(data = Data, aes(x = V1,y = V2)) +  geom_bar(stat="identity") +
  labs(x = 'Fold change',y = 'Frequency',title = paste('Module',i)) + 
  theme(plot.title=element_text(vjust=0.5,hjust=0.5))
}
plot_grid(plotlist = g, nrow = 2)

write.csv(gene_module,'~/Downloads/genes.csv')
TF_module <- net$colors[TF$gene[TF$cluster==2]]

i=1;{
enrichR <- read.table(paste('~/Downloads/enrichR/module',i,'.txt',sep=''),sep ="\t",header = T)
enrichR$TF <- str_split_fixed(enrichR$Term,"_",2)[,1]
#EnrichTF<-enrichR[1:30,]
EnrichTF <- enrichR[enrichR$TF%in%intersect(enrichR$TF,toupper(names(TF_module)[TF_module==i])),]
EnrichTF$geneCount <-as.numeric(str_split_fixed(EnrichTF$Overlap,"/",2)[,1])
EnrichTF$geneTotal <-as.numeric(str_split_fixed(EnrichTF$Overlap,"/",2)[,2])
#which(rownames(EnrichTF)=="50")
#EnrichTF<-EnrichTF[-13,]

EnrichTF$TF <- factor(EnrichTF$TF, levels=rev(unique(EnrichTF$TF)))

ggplot(data = EnrichTF, aes(x = EnrichTF$Combined.Score,y = EnrichTF$TF),) +
  geom_point(aes(size = EnrichTF$geneCount,
                 color = -1*log10(EnrichTF$P.value))) +
  scale_color_gradient(low="green",high="red") +
  labs(color = expression(-log[10](P.value)),size = "Gene Count",
       x = "Combined Score",y="Transcriptional factor",
       title=paste('TF enrichment of module', i ,'DEG')) +
  theme(plot.title=element_text(size=20,vjust=0.5,hjust=0.5),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        legend.position="right",
        axis.text=element_text(size=15),
        axis.title=element_text(size=18),)
}

splenicB <- read.xlsx('~/Dropbox/77-Tcf1 in naive CD8/introduction/DEG/splenic B cells.xlsx')



#DEG calculation
HdacD01.DEG <- FindAllMarkers(object = HdacD01)
write.csv(HdacD01.markers,"~/R/scRNAseq_raw/All DEG.csv")
markers3 <- FindMarkers(object=HdacD01,ident.1=3)
TFmarker1<-markers1[intersect(TransFactor$Symbol,rownames(markers1)),]

markers_D0<-FindMarkers(object=HdacD01,ident.1=HdacD0.data@Dimnames[[2]],
                     ident.2=setdiff(HdacD0_Dko.data@Dimnames[[2]],c("CTCACACAGAACTGTA-1")))
markers_D1<-FindMarkers(object=HdacD01,ident.1=HdacD1.data@Dimnames[[2]],
                          ident.2=HdacD1_Dko.data@Dimnames[[2]])
markers_Ctrl<-FindMarkers(object=HdacD01,ident.1=HdacD0.data@Dimnames[[2]],
                        ident.2=setdiff(HdacD1.data@Dimnames[[2]],c("CTCACACAGAACTGTA-1")))
markers_Dko<-FindMarkers(object=HdacD01,ident.1=HdacD0.data@Dimnames[[2]],
                        ident.2=HdacD1_Dko.data@Dimnames[[2]])

#==========================Ctrl Div 1 subset DEG analysis

HdacD1 <- CreateSeuratObject(counts = as.matrix(HdacD1.data)[AllMarkers$gene[AllMarkers$cluster==2],], 
                              min.cells = 3, min.features = 200,
                              project = "10X_HdacD0")
HdacD1 <- NormalizeData(object = HdacD1, 
                         normalization.method = "LogNormalize",
                         scale.factor = 10000)
HdacD1 <- FindVariableFeatures(object = HdacD1, 
                                mean.function = ExpMean, 
                                dispersion.function = LogVMR, 
                                x.low.cutoff = 0, 
                                x.high.cutoff = 4, 
                                y.cutoff = 0.5)
VariableFeaturePlot(object = HdacD1)
all.genes<-rownames(HdacD1)
HdacD1 <- ScaleData(object = HdacD1, features = all.genes)
HdacD1 <- RunPCA(object = HdacD1, features = VariableFeatures(object = HdacD1))
HdacD1 <- FindNeighbors(object = HdacD1, dims = 1:10)
HdacD1<-FindClusters(object = HdacD1,resolution = 2)
HdacD1<-RunTSNE(object = HdacD1, dims.use = 1:10)
DimPlot(HdacD1)

#Find differentially expressed genes
markersD1 <- FindAllMarkers(HdacD1,only.pos=F)
FeaturePlot(HdacD1,markersD1$gene[markersD1$cluster==2][1:9])
write.csv(markersD1,'~/Downloads/markers of each subcluster.csv')
markersD1$gene[markersD1$cluster==2 & markersD1$avg_logFC>0]


ExprMtxD1 = as.data.frame(t(as.matrix(HdacD1@assays$RNA@data[AllMarkers$gene[AllMarkers$cluster==2],])))
ExprMtxD1$ident <- as.numeric(HdacD1@active.ident)
ExpMean <- aggregate(ExprMtxD1,list(ExprMtxD1$ident),mean)[,c(-1,-1044)]

geneCorD1 <- cor(ExprMean)
geneCorD1 <- geneCorD1[which(is.na(geneCorD1[,1])==F),which(is.na(geneCorD1[,1])==F)]
geneCorD1[which(geneCorD1==1)]<-0
hcluter=hclust(dist(geneCorD1))
callback = function(hcluter, geneCorD1){as.hclust(reorder(hcluter,dist(geneCorD1),method = "OLO"))}
pheatmap(geneCorD1,clustering_callback=callback,cutree_rows = 6, cutree_cols = 6, show_colnames = F,show_rownames = F)
hist(geneCorD1)


#==========================correlation with Myc
corMyc = 0;
for (i in 1:12825){
corMyc[i] <- cor(as.numeric(HdacD01@assays$RNA@data["Myc",c(HdacD0.data@Dimnames[[2]],HdacD1.data@Dimnames[[2]])]),
    as.numeric(HdacD01@assays$RNA@data[i,c(HdacD0.data@Dimnames[[2]],HdacD1.data@Dimnames[[2]])]))
}
CorMyc <- as.data.frame(cbind(rownames(HdacD01@assays$RNA@data),corMyc))
CorMyc$corMyc <- corMyc; CorMyc$V1 <- as.character(CorMyc$V1)
CorMyc <-na.omit(CorMyc)
CorMyc<-CorMyc[order(CorMyc$corMyc),]



ExprMtx = as.data.frame(t(as.matrix(HdacD01@assays$RNA@data[AllMarkers$gene[AllMarkers$cluster==2],
                                                            HdacD1.data@Dimnames[[2]]])))

write.csv(geneCor,file='~/Downloads/corr.csv')

geneCor <- cor(ExprMtx)
geneCor <- geneCor[which(is.na(geneCor[,1])==F),which(is.na(geneCor[,1])==F)]
geneCor[which(geneCor==1)]<-0
Heatmap <- pheatmap(geneCor,show_colnames = F,show_rownames = F)#,kmeans_k = 4,
                    #cluster_cols=F,cluster_rows=F)
hcluter=hclust(dist(geneCor))
callback = function(hcluter, geneCor){as.hclust(reorder(hcluter,dist(geneCor),method = "OLO"))}
geneCor[which(geneCor>0.5)]<-0
Heatmap1 <- pheatmap(geneCor,clustering_callback=callback,cutree_rows = 5, cutree_cols = 5, show_colnames = F,show_rownames = F)



gene_cluster_name[(1+107+377+464):(107+377+464+95)] = names(Heatmap$kmeans$cluster)[Heatmap$kmeans$cluster==4]

#gene_cluster_name[(96+221+422):(95+221+422+305)] = names(sort(apply(gene_cluster_cor[names(Heatmap$kmeans$cluster)[Heatmap$kmeans$cluster==1],],1,mean),decreasing = T))

data("eurodist")
d <- as.dist(eurodist)
hc <- hclust(dist(eurodist))

## plot original dendrogram and the reordered dendrograms  
plot(hc)  
plot(reorder(hc, d, method = "GW"))  
plot(reorder(hc, d, method = "OLO"))

intersect(names(Heatmap$kmeans$cluster)[Heatmap$kmeans$cluster==1],AllMarkers$gene[AllMarkers$cluster==2][1:107])
gene_cluster_cor = cor(ExprMtx[,gene_cluster_name[(1+107+377+464):(107+377+464+95)]])
gene_cluster_cor[which(gene_cluster_cor==1)]<-0
hmap = pheatmap(gene_cluster_cor,show_colnames = F,show_rownames = F,cluster_rows=T,cluster_cols=T)
                #gaps_row=c(107,107+377,107+377+464),gaps_col=c(107,107+377,107+377+464))


write.csv(Seq,"gene.csv")
which(Seq%in%c("E2f4","Ssrp1","Atf1","Mta2","Rbck1","Ubtf","Yy1","Ybx1"))

write.csv(gene_mean,'~/R/scRNAseq_raw/marker_mean_Expr.csv')

#scatter plot 
gene_cluster=Heatmap1$tree_row$labels[Heatmap1$tree_row$order]#names(Heatmap$kmeans$cluster)[Heatmap$kmeans$cluster==3]
Expr_Mtx = as.data.frame(as.matrix(HdacD01@assays$RNA@data[gene_cluster,]))
expr_mean <- as.data.frame(matrix(NA,nrow=nrow(Expr_Mtx),ncol=4),row.names=rownames(Expr_Mtx))
names(expr_mean) = c('Ctrl0','Dko0','Ctrl1','Dko1')
expr_mean$Ctrl0 <- apply(Expr_Mtx[,HdacD0.data@Dimnames[[2]]],1,mean)
expr_mean$Ctrl1 <- apply(Expr_Mtx[,HdacD1.data@Dimnames[[2]]],1,mean)
expr_mean$Dko0 <- apply(Expr_Mtx[,HdacD0_Dko.data@Dimnames[[2]]],1,mean)
expr_mean$Dko1 <- apply(Expr_Mtx[,HdacD1_Dko.data@Dimnames[[2]]],1,mean)
expr_mean$genes <- rownames(expr_mean)

ggplot(data=NULL,aes(x = 1:length(expr_mean$Ctrl1),y = expr_mean$Ctrl1))+geom_line()+ 
  labs(x = "gene",y="mean expression",title="gene mean expression in cluster 2")+
  theme(plot.title=element_text(vjust=0.5,hjust=0.5))


ggplot(data=expr_mean,aes(x = Dko1,y = Ctrl1)) + 
  geom_point(color='blue') + labs(x = expression(Dko1),y=expression(Ctrl1),title="Scatter plot of group 4 genes") +
  geom_smooth(method='lm') + geom_text_repel(data=subset(expr_mean,Ctrl1-Dko1<0),aes(label=genes,hjust=2)) + 
  #geom_text_repel(data=subset(expr_mean,expr_mean$genes %in% TF$gene),aes(label=genes,hjust=-2,colour='red'),show.legend=F) + 
  geom_line(aes(x = seq(0,max(expr_mean$Ctrl1),length.out=nrow(expr_mean)),y=seq(0,max(expr_mean$Ctrl1),length.out=nrow(expr_mean))))+
  theme(plot.title=element_text(vjust=0.5,hjust=0.5))
intersect(gene_cluster,TF$gene)

ggplot(data=expr_mean,aes(x = Dko0,y = Dko1)) + 
  geom_point(color='blue') + labs(x = "Dko0",y="Dko1",title="Scatter plot of group 4 genes") +
  geom_smooth(method='lm') + geom_text_repel(data=subset(expr_mean,abs(Dko1-Dko0)>0.3),aes(label=genes,vjust=-4)) + 
  #geom_text_repel(data=subset(expr_mean,expr_mean$genes %in% TF$gene),aes(label=genes,hjust=-2,colour='red'),show.legend=F) + 
  geom_line(aes(x = seq(0,max(expr_mean$Dko1),length.out=nrow(expr_mean)),y=seq(0,max(expr_mean$Dko1),length.out=nrow(expr_mean))))+
  theme(plot.title=element_text(vjust=0.5,hjust=0.5))

ggplot(data=expr_mean,aes(x = Ctrl0,y = Ctrl1)) + 
  geom_point(color='blue') + labs(x = "Ctrl0",y="Ctrl1",title="Scatter plot of group 4 genes") +
  geom_smooth(method='lm') + geom_text_repel(data=subset(expr_mean,abs(Ctrl1-Ctrl0)>0.7),aes(label=genes,vjust=-3)) + 
  #geom_text_repel(data=subset(expr_mean,expr_mean$genes %in% TF$gene),aes(label=genes,hjust=-2,colour='red'),show.legend=F) + 
  geom_line(aes(x = seq(0,max(expr_mean$Ctrl1),length.out=nrow(expr_mean)),y=seq(0,max(expr_mean$Ctrl1),length.out=nrow(expr_mean))))+
  theme(plot.title=element_text(vjust=0.5,hjust=0.5))

ggplot(data=expr_mean,aes(x = Ctrl0,y = Dko0)) + 
  geom_point(color='blue') + labs(x = "Ctrl0",y="Dko0",title="Scatter plot of group 4 genes") +
  geom_smooth(method='lm') + geom_text_repel(data=subset(expr_mean,abs(Dko0-Ctrl0)>0.2),aes(label=genes,hjust=-1)) + 
  #geom_text_repel(data=subset(expr_mean,expr_mean$genes %in% TF$gene),aes(label=genes,hjust=-2,colour='red'),show.legend=F) + 
  geom_line(aes(x = seq(0,max(expr_mean$Dko0),length.out=nrow(expr_mean)),y=seq(0,max(expr_mean$Dko0),length.out=nrow(expr_mean))))+
  theme(plot.title=element_text(vjust=0.5,hjust=0.5))



intersect(gene_cluster_cor,TF$gene)

gene_cluster=names(Heatmap$kmeans$cluster)[Heatmap$kmeans$cluster==4]
write.csv(gene_cluster,"~/Downloads/gene.csv")

a = read.csv('~/Downloads/BigWig/macs2/venn/chipseq_no_cutrun.txt')
b = read.csv('~/Downloads/BigWig/macs2/venn/chipseq_cutrun.txt')

ExprMtx$cell <- as.factor(ifelse(rownames(ExprMtx)%in%HdacD0.data@Dimnames[[2]],"Div0","Div1"))
ZeroCount <- which(ExprMtx[ExprMtx$Myc==0,rev(CorMyc$V1)[rev(CorMyc$corMyc)>0.15]]==0,arr.ind = T)
x <- as.data.frame(sort(table(ZeroCount[,2]),decreasing = T))
#ggplot(data=NULL,aes(x = 1:(dim(CorMyc)[1]-1),y = CorMyc$corMyc[1:(dim(CorMyc)[1]-1)])) + 
g = 1:20; gene = rev(CorMyc$V1)[2:21]#rev(CorMyc$V1)[as.numeric(as.character(x$Var1[g]))]; 
ggplot(data=ExprMtx,aes(x = ExprMtx$Myc,y = ExprMtx[,gene])) + 
  geom_point(aes(color = ExprMtx$cell)) + 
  scale_color_manual(values=c("Div0"="red","Div1"="blue"))+
  labs(x="Myc",y= gene,color="Cell type",
       title = paste(601-x$Freq[g],"    cor",round(cor(ExprMtx$Myc,ExprMtx[,gene]),2))) +
  #geom_hline(yintercept= 0,linetype=4) + 
  theme(plot.title=element_text(size=20),#,vjust=0.5,hjust=0.5),
      legend.title=element_text(size=15),
      legend.text=element_text(size=15),
      legend.position="right",
      axis.text=element_text(size=15),
      axis.title=element_text(size=18))



#for (i in as.character(CorMyc$V1[1:9])){
FeatureScatter(HdacD01,feature1 = "Myc", feature2 = "Mybl1",#as.character(CorMyc$V1)[1],slot = "data",
               cells=c(HdacD0.data@Dimnames[[2]],HdacD1.data@Dimnames[[2]]))
FeatureScatter(HdacD01,feature1 = "Myc", feature2 = as.character(rev(CorMyc$V1))[2],slot = "data",
               cells=c(HdacD0.data@Dimnames[[2]],HdacD1.data@Dimnames[[2]]))

#calculating markers
write.csv(intersect(TransFactor$Symbol,rownames(HdacDEG)),"~/R/scRNAseq_raw/TFgene.csv")

HdacDEG <- markers_Ctrl[intersect(rownames(markers_Ctrl),rownames(markers_Dko)),]
NonchangeDEG <- markers_D0[intersect(rownames(markers_D0),rownames(markers_D1)),]
write.csv(NonchangeDEG,"~/R/scRNAseq_raw/NonchangeDEG from D0 to D1.csv")


top10 <- HdacD01.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(HdacD01, features = top10$gene) + NoLegend()

HdacD01_subset <- subset(x=HdacD01, cells=HdacD01.data@Dimnames[[2]][HdacD01@active.ident==3])

FeaturePlot(#HdacD01_subset,AllMarkers$gene[AllMarkers$cluster==3][1:9],
            #HdacD01,HdacD01.markers$gene[HdacD01.markers$cluster==3][1:9],
            #HdacD01_subset,c("mt-Nd4","mt-Cytb","mt-Nd5","mt-Atp8","Chpt1","Rrm2b","Gm29666","B930036N10Rik","Aldh1a3"),
            HdacD01,HdacD01.markers$gene[HdacD01.markers$cluster==3][1:9],
            #c("Myc", "Egr1","Egr2","Egr3","Jun","Fosb"),
            #c("mt-Co1","Actb","mt-Nd4",rownames(markers3)[1:3]),
            #rownames(markers2)[c(4,6:10)],
            #rownames(markers1)[c(1:3,6,10,12)],
            reduction="umap",ncol=3,cols = c("green", "red"),label=T)

write.csv(t(TF$gene[TF$cluster==2]),file = "~/R/scRNAseq_raw/TFgene.csv")

FeaturePlot(HdacD01,top10$gene[top10$cluster==3][1:9],reduction="tsne",
            cols = c("green", "red"),ncol=3,label=T,label.size = 5)

RidgePlot(HdacD01,HdacD01.markers$gene[HdacD01.markers$cluster==2][1:9],
            #HdacD01,rownames(TFmarker1)[TFmarker1$avg_logFC<0],
            #HdacD01_subset,c("mt-Nd4","mt-Cytb","mt-Nd5","mt-Atp8","Chpt1","Rrm2b","Gm29666","B930036N10Rik","Aldh1a3"),
            #HdacD01_subset,AllMarkers$gene[AllMarkers$cluster==3][1:9],
            #HdacD01, features = TF$gene[TF$cluster==3],
            #c("Myc", "Egr1","Egr2","Egr3","Jun","Fosb"),
            #rownames(markers)[1:4],
          ncol = 3)

TransFactor <- read.table("~/Downloads/genomicData/Mus_musculus_TF.txt",sep ="\t",header = T)
TF <- AllMarkers[AllMarkers$gene%in%intersect(TransFactor$Symbol,AllMarkers$gene),]

housekeeper = c("Rps18","Gapdh","Rplp0","B2m")
mean(HdacD0@assays$RNA@counts["Gapdh",])
sd(Hdac@assays$RNA@counts["Rps18",])

intersect(read.csv("~/R/scRNAseq_raw/D0_Dko_vs_D0_Ctrl.csv")[,"X"],TransFactor$Symbol)

enrichR <- read.table("~/Downloads/enrichR/module1.txt",sep ="\t",header = T)
enrichR$TF <- str_split_fixed(enrichR$Term," ",2)[,1]
EnrichTF<-enrichR[1:30,]
EnrichTF <- enrichR[enrichR$TF%in%intersect(enrichR$TF,toupper(TF$gene)),]
EnrichTF$geneCount <-as.numeric(str_split_fixed(EnrichTF$Overlap,"/",2)[,1])
EnrichTF$geneTotal <-as.numeric(str_split_fixed(EnrichTF$Overlap,"/",2)[,2])
#which(rownames(EnrichTF)=="50")
#EnrichTF<-EnrichTF[-13,]

EnrichTF$TF <- factor(EnrichTF$TF, levels=rev(unique(EnrichTF$TF)))

ggplot(data = EnrichTF, aes(x = EnrichTF$Combined.Score,y = EnrichTF$TF),) +
  geom_point(aes(size = EnrichTF$geneCount,
                 color = -1*log10(EnrichTF$P.value))) +
  scale_color_gradient(low="green",high="red") +
  labs(color = expression(-log[10](P.value)),size = "Gene Count",
       x = "Combined Score",y="Epigenomic site",
       title="Epigenomic enrichment of cluster 1 DEG") +
  theme(plot.title=element_text(size=20,vjust=0.5,hjust=0.5),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        legend.position="right",
        axis.text=element_text(size=15),
        axis.title=element_text(size=18),)


library("sva")
library(bladderbatch)
data(bladderdata)
dat <- bladderEset[1:50,]

pheno = pData(dat)
edata = exprs(dat)
batch = pheno$batch
mod = model.matrix(~as.factor(cancer), data=pheno)

ExprAll = as.matrix(Integrated@assays$RNA@counts)
batch = as.numeric(str_split_fixed(colnames(ExprAll),"_",2)[,2])

# parametric adjustment
combat_ExprAll = ComBat(dat=ExprAll, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)

# non-parametric adjustment, mean-only version
combat_edata2 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=FALSE, mean.only=TRUE)

# reference-batch version, with covariates
combat_edata3 = ComBat(dat=edata, batch=batch, mod=mod, par.prior=TRUE, ref.batch=3)
