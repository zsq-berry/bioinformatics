library(BiocParallel)
register(MulticoreParam(8))
library(chromVAR)
library(motifmatchr)
library(BSgenome.Mmusculus.UCSC.mm9)
library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomicRanges)
library(SummarizedExperiment)
library(ggplot2)
library(ggrepel)
library(stringr)
library(pheatmap)
library(devtools)
library(TFBSTools)
### Using example counts from package ------------------------------------------

install("~/Downloads/chromVARmotifs-master/")
load("~/Downloads/chromVARmotifs-master/data/mouse_pwms_v1.rda")
load("~/Downloads/chromVARmotifs-master/data/mouse_pwms_v2.rda")
load("~/Downloads/chromVARmotifs-master/data/human_pwms_v1.rda")
coverage <- as.matrix(coverage)

{
fragment_counts <- SummarizedExperiment(assays = list(counts = coverage),
                                        rowRanges =GRanges(merge))
fragment_counts <- addGCBias(fragment_counts, 
                            genome = BSgenome.Mmusculus.UCSC.mm9)
# counts_filtered <- filterSamples(fragment_counts, min_depth = 1500,
#                                  min_in_peaks = 0.15)
# counts_filtered <- filterPeaks(fragment_counts)
motif_ix <- matchMotifs(mouse_pwms_v1, fragment_counts,
                        genome = BSgenome.Mmusculus.UCSC.mm9)

# computing deviations
dev <- computeDeviations(object = fragment_counts, 
                         annotations = motif_ix)

variability <- computeVariability(dev)

plotVariability(variability, use_plotly = FALSE, n = 10) 

tsne_results <- deviationsTsne(dev, threshold = 1.5, perplexity = 10)

tsne_plots <- plotDeviationsTsne(dev, tsne_results,
                                 annotation_name = "Tcf7",
                                 #sample_column = "Cell_Type",
                                 shiny = FALSE)

tsne <- tsne_plots[[1]]$data

# ggplot(data = tsne,aes(x=x, y=y, label=text, colour = color)) +
#   geom_point() + geom_text_repel(aes(label=text), size = 4, hjust=0, vjust=-1, angle=0) +
#   scale_color_gradient(low = "green", high = "red") +
#   labs(color = 'Tcf7',x = "tsne1",y="tsne2")

# order(variability$variability,decreasing = T)
# 
# variability$name[order(variability$variability,decreasing = T)[1:10]]
# which(variability$name=='Ctcf')

rownames(dev@assays@data$z) <- as.character(motif_ix$name)
deviation <- dev@assays@data$z[order(variability$variability,decreasing = T),rename[1:10]]
#rownames(deviation) <- str_split_fixed(rownames(deviation),'_',2)[,2]
pheatmap(deviation[1:50,],cluster_cols = F,cluster_rows = F,show_rownames = T,angle_col='315',
         main= paste('TF deviation z score in cluster', 'all'))
}
cbind(deviation['Ctcf',],1:10)

motif_in_cluster <- data.frame(matrix(data=NA,nrow=20,ncol=14))
#motif_match_num <- matrix(data=NA,nrow=50,ncol=7)
for (i in 1) {
  union_idx <- as.numeric(rownames(diff_peak_kmeans)[diff_peak_kmeans$cluster_idx==i])
  frag_counts <- SummarizedExperiment(assays = list(counts = coverage[union_idx,]),
                                          rowRanges =GRanges(merge[union_idx,]))
  frag_counts <- addGCBias(frag_counts, 
                               genome = BSgenome.Mmusculus.UCSC.mm9)
  motif_ix <- matchMotifs(mouse_pwms_v1, frag_counts,
                          genome = BSgenome.Mmusculus.UCSC.mm9)
  
  match_num <- apply(motif_ix@assays@data$motifMatches,2,sum)
  
  motif_in_cluster[,2*i-1] <- motif_ix@colData[order(match_num,decreasing = T)[1:20],]
  
  motif_in_cluster[,2*i] <- sort(match_num,decreasing = T)[1:20]/cluster_num[i]
  #motif_match_num[,i] <- sort(match_num,decreasing = T)[1:20]
}

as.character(motif_ix@colData[order(match_num,decreasing = T)[1:20],])


for(i in 1:7){
  pheatmap(deviation[motif_in_cluster[,i],],cluster_cols = F,cluster_rows = T,show_rownames = T,angle_col='315',
           main= paste('TF deviation z score in cluster', i))
}

motif_name <- as.data.frame(motif_ix@colData$name)


hits <-suppressWarnings(findOverlaps(GRanges(diff_peak_kmeans),GRanges(mm9_promoter)))
cbind(1:7,table(diff_peak_kmeans[hits@from,"cluster_idx"]))


ctcf_motif <- ctcf_motif <- matchMotifs(mouse_pwms_v1, GRanges(CTCF_B_Blood

ctcf_motif_peak <- subset(diff_peak_kmeans,cluster_idx==1)[ctcf_motif,]

hits <- findOverlaps(GRanges(tcf1),GRanges(ctcf_motif_peak))
hits1 <- findOverlaps(GRanges(ctcf_motif_peak[hits@to,]),GRanges(mm9_deg_50kb))
cbind(1:8,table(mm9_deg[hits1@to,]$cluster))
mm9_deg[hits1@to,]$gene_name[mm9_deg[hits1@to,]$cluster==2]

stat5_wt_24h_summit <- stat5_wt_24h
stat5_wt_24h_summit$start <- stat5_wt_24h$start+stat5_wt_24h$V10-150
stat5_wt_24h_summit$end <- stat5_wt_24h$start+stat5_wt_24h$V10+150
stat5_wt_24h_summit <- stat5_wt_24h_summit[-which(stat5_wt_24h_summit$chr=='chrX'),]

stat5_dko_24h_summit <- stat5_dko_24h
stat5_dko_24h_summit$start <- stat5_dko_24h$start+stat5_dko_24h$V10-150
stat5_dko_24h_summit$end <- stat5_dko_24h$start+stat5_dko_24h$V10+150
stat5_dko_24h_summit <- stat5_dko_24h_summit[-which(stat5_dko_24h_summit$chr=='chrX'),]

motif_idx <- matchMotifs(mouse_pwms_v1, GRanges(tcf1[order(tcf1$V5)[1:1000],]),
            genome = BSgenome.Mmusculus.UCSC.mm9)

match_num <- apply(motif_idx@assays@data$motifMatches,2,sum)

as.data.frame(cbind(as.character(motif_idx@colData[order(match_num,decreasing = T)[1:100],]),
                    sort(match_num,decreasing = T)[1:100]))


match_num[which(motif_ix@colData$name=='Tcf7')]
