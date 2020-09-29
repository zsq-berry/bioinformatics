Tcf1_liver <- read.table('/media/shaoqizhu/easystore/pioneer_factor/HepG2_liver_Epithelium_Tcf1_ENCSR444LIN_2.bed')
Tcf1_boneMarrow <- read.table('/media/shaoqizhu/easystore/pioneer_factor/K562_boneMarrow_Tcf1_ENCSR863KUB_2.bed')
Ctcf_liver <- read.table('/media/shaoqizhu/easystore/pioneer_factor/HepG2_liver_Epithemium_Ctcf_GSM749683.bed')
Ctcf_boneMarrow <- read.table('/media/shaoqizhu/easystore/pioneer_factor/K562_boneMarrow_Ctcf_GSM624080.bed')
colnames(Tcf1_liver)[1:3] <- c('chr','start','end')
colnames(Tcf1_boneMarrow)[1:3] <- c('chr','start','end')
colnames(Ctcf_liver)[1:3] <- c('chr','start','end')
colnames(Ctcf_boneMarrow)[1:3] <- c('chr','start','end')

hits1 <- findOverlaps(GRanges(Tcf1_liver),GRanges(Ctcf_liver))
hits2 <- findOverlaps(GRanges(Tcf1_boneMarrow),GRanges(Ctcf_boneMarrow))
hits1

BATF_BLymph <- read.table('/media/shaoqizhu/easystore/pioneer_factor/Homo_BLymph_blood_BATF_GSM803538.bed')
CTCF_BLymph <- read.table('/media/shaoqizhu/easystore/pioneer_factor/Homo_BLymph_blood_Ctcf_GSM1233888.bed')
colnames(BATF_BLymph)[1:3] <- c('chr','start','end')
colnames(CTCF_BLymph)[1:3] <- c('chr','start','end')
overlap <- findOverlaps(GRanges(BATF_BLymph),GRanges(CTCF_BLymph))
draw.pairwise.venn(nrow(BATF_BLymph),nrow(CTCF_BLymph),length(overlap@from))


BRN2_Neur <- read.table('/media/shaoqizhu/easystore/pioneer_factor/Homo_NeuralProgen_BRN2_GSM1934425.bed')
CTCF_Neur <- read.table('/media/shaoqizhu/easystore/pioneer_factor/Homo_NeuralProgen_CTCF_GSM2259909.bed')
colnames(BRN2_Neur)[1:3] <- c('chr','start','end')
colnames(CTCF_Neur)[1:3] <- c('chr','start','end')
findOverlaps(GRanges(BRN2_Neur),GRanges(CTCF_Neur))

REST_Neur <- read.table('/media/shaoqizhu/easystore/pioneer_factor/Neural_Progen_REST_GSM671100.bed')
CTCF_Neur <- read.table('/media/shaoqizhu/easystore/pioneer_factor/Neural_Progen_Ctcf_GSM671100.bed')
colnames(REST_Neur)[1:3] <- c('chr','start','end')
colnames(CTCF_Neur)[1:3] <- c('chr','start','end')
findOverlaps(GRanges(REST_Neur),GRanges(CTCF_Neur))

CEBPA_MP <- read.table('/media/shaoqizhu/easystore/pioneer_factor/Mus_Macrophage_BM_CEBPA_GSM2845729.bed')
CTCF_MP <- read.table('/media/shaoqizhu/easystore/pioneer_factor/Mus_Macrophage_BM_CTCF_GSM918726.bed')
colnames(CEBPA_MP)[1:3] <- c('chr','start','end')
colnames(CTCF_MP)[1:3] <- c('chr','start','end')
findOverlaps(GRanges(CEBPA_MP),GRanges(CTCF_MP))

FoxA1_Breast <- read.table('/media/shaoqizhu/easystore/pioneer_factor/Homo_Breast_FoxA1_GSM1858629.bed')
CTCF_Breast <- read.table('/media/shaoqizhu/easystore/pioneer_factor/Homo_Breast_Ctcf_GSM3357135.bed')
colnames(FoxA1_Breast)[1:3] <- c('chr','start','end')
colnames(CTCF_Breast)[1:3] <- c('chr','start','end')
findOverlaps(GRanges(FoxA1_Breast),GRanges(CTCF_Breast))

Gata2_Prostate <- read.table('/media/shaoqizhu/easystore/pioneer_factor/Homo_Prostate_GATA2_LNCaP_GSM1691145.bed')
CTCF_Prostate <- read.table('/media/shaoqizhu/easystore/pioneer_factor/Homo_Prostate_CTCF_LNCaP_GSM2827203.bed')
colnames(Gata2_Prostate)[1:3] <- c('chr','start','end')
colnames(CTCF_Prostate)[1:3] <- c('chr','start','end')
findOverlaps(GRanges(Gata2_Prostate),GRanges(CTCF_Prostate))

Pax7_corticotroph <- read.table('/media/shaoqizhu/easystore/pioneer_factor/Mus_Corticotroph_Pax7_GSM2445273.bed')
CTCF_corticotroph <- read.table('/media/shaoqizhu/easystore/pioneer_factor/Mus_corticotroph_Ctcf_GSM2684854.bed')
colnames(Pax7_corticotroph)[1:3] <- c('chr','start','end')
colnames(CTCF_corticotroph)[1:3] <- c('chr','start','end')
findOverlaps(GRanges(Pax7_corticotroph),GRanges(CTCF_corticotroph))

TP63_Keratinocyte <- read.table('/media/shaoqizhu/easystore/pioneer_factor/Homo_Keratinocyte_TP63_GSM1366694.bed')
CTCF_Keratinocyte <- read.table('/media/shaoqizhu/easystore/pioneer_factor/Homo_Keratinocyte_CTCF_GSM733636.bed')
colnames(TP63_Keratinocyte)[1:3] <- c('chr','start','end')
colnames(CTCF_Keratinocyte)[1:3] <- c('chr','start','end')
findOverlaps(GRanges(TP63_Keratinocyte),GRanges(CTCF_Keratinocyte))

#Cell line GM12878
B_Blood_TF <- read.table('/media/shaoqizhu/easystore/pioneer_factor/B_Blood_Pax5_GSM803391.bed')
B_Blood_CTCF <- read.table('/media/shaoqizhu/easystore/pioneer_factor/B_Blood_Ctcf_GSM935611.bed')
colnames(B_Blood_TF)[1:3] <- c('chr','start','end')
colnames(B_Blood_CTCF)[1:3] <- c('chr','start','end')
overlap <- findOverlaps(GRanges(B_Blood_TF),GRanges(B_Blood_CTCF))
draw.pairwise.venn(nrow(B_Blood_TF),nrow(B_Blood_CTCF),length(overlap@from))

#Cell line IMR90
TF_Lung_fibro <- read.table('/media/shaoqizhu/easystore/pioneer_factor/Lung_fibroblast_RELA_GSM1055810.bed')
CTCF_Lung_fibro <- read.table('/media/shaoqizhu/easystore/pioneer_factor/Lung_fibroblast_CTCF_ GSM935404.bed')
colnames(TF_Lung_fibro)[1:3] <- c('chr','start','end')
colnames(CTCF_Lung_fibro)[1:3] <- c('chr','start','end')
overlap <- findOverlaps(GRanges(TF_Lung_fibro),GRanges(CTCF_Lung_fibro))
draw.pairwise.venn(nrow(TF_Lung_fibro),nrow(CTCF_Lung_fibro),length(overlap@from))

#Cell line MCF-7
Breast_epithelium_TF <- read.table('/media/shaoqizhu/easystore/pioneer_factor/Breast_epithelium_TFAP2A_GSM588928.bed')
Breast_epithelium_CTCF <- read.table('/media/shaoqizhu/easystore/pioneer_factor/Breast_epithelium_CTCF_GSM614615.bed')
colnames(Breast_epithelium_TF)[1:3] <- c('chr','start','end')
colnames(Breast_epithelium_CTCF)[1:3] <- c('chr','start','end')
overlap <- findOverlaps(GRanges(Breast_epithelium_TF),GRanges(Breast_epithelium_CTCF))
draw.pairwise.venn(nrow(Breast_epithelium_TF),nrow(Breast_epithelium_CTCF),length(overlap@from))

#Cell line HepG2
TF_Liver_epithelium <- read.table('/media/shaoqizhu/easystore/pioneer_factor/Liver_epithelium_NFIC_GSM1010741.bed')
CTCF_Liver_epithelium <- read.table('/media/shaoqizhu/easystore/pioneer_factor/Liver_epithelium_CTCF_GSM749683.bed')
colnames(TF_Liver_epithelium)[1:3] <- c('chr','start','end')
colnames(CTCF_Liver_epithelium)[1:3] <- c('chr','start','end')
overlap <- findOverlaps(GRanges(TF_Liver_epithelium),GRanges(CTCF_Liver_epithelium))
draw.pairwise.venn(nrow(TF_Liver_epithelium),nrow(CTCF_Liver_epithelium),length(overlap@from))

#Cell line K562
Bone_marrow_TF <- read.table('/media/shaoqizhu/easystore/pioneer_factor/Bone_marrow_NFYA_GSM648587.bed')
Bone_marrow_CTCF <- read.table('/media/shaoqizhu/easystore/pioneer_factor/Bone_marrow_CTCF_GSM646432.bed')
colnames(Bone_marrow_TF)[1:3] <- c('chr','start','end')
colnames(Bone_marrow_CTCF)[1:3] <- c('chr','start','end')
overlap <- findOverlaps(GRanges(Bone_marrow_TF),GRanges(Bone_marrow_CTCF))
draw.pairwise.venn(nrow(Bone_marrow_TF),nrow(Bone_marrow_CTCF),length(overlap@from))

#Cell line SK-N-SH
TF_Neuroblastoma <- read.table('/media/shaoqizhu/easystore/pioneer_factor/Neuroblastoma_RARA_GSM1693101.bed')
CTCF_Neuroblastoma <- read.table('/media/shaoqizhu/easystore/pioneer_factor/Neuroblastoma_CTCF_GSM3498398.bed')
colnames(TF_Neuroblastoma)[1:3] <- c('chr','start','end')
colnames(CTCF_Neuroblastoma)[1:3] <- c('chr','start','end')
overlap <- findOverlaps(GRanges(TF_Neuroblastoma),GRanges(CTCF_Neuroblastoma))
draw.pairwise.venn(nrow(TF_Neuroblastoma),nrow(CTCF_Neuroblastoma),length(overlap@from))

Gata3_T_Blood <- read.table('/media/shaoqizhu/easystore/pioneer_factor/Mus_T_Blood_Gata3_GSM523229.bed')
CTCF_T_Blood <- read.table('/media/shaoqizhu/easystore/pioneer_factor/Mus_T_Blood_Ctcf_GSM1023416.bed')
colnames(Gata3_T_Blood)[1:3] <- c('chr','start','end')
colnames(CTCF_T_Blood)[1:3] <- c('chr','start','end')
findOverlaps(GRanges(Gata3_T_Blood),GRanges(CTCF_T_Blood))

Myod1_Muscle <- read.table('/media/shaoqizhu/easystore/pioneer_factor/Mus_Muscle_Myod1_GSM539539.bed')
CTCF_Muscle <- read.table('/media/shaoqizhu/easystore/pioneer_factor/Mus_Muscle_Ctcf_GSM915188.bed')
colnames(Myod1_Muscle)[1:3] <- c('chr','start','end')
colnames(CTCF_Muscle)[1:3] <- c('chr','start','end')
findOverlaps(GRanges(Myod1_Muscle),GRanges(CTCF_Muscle))


hits1 <- findOverlaps(GRanges(CTCF_Bone_marrow),GRanges(CTCF_B_Blood))
hits2 <- findOverlaps(GRanges(CTCF_B_Blood),GRanges(CTCF_Liver_epithelium))
hits3 <- findOverlaps(GRanges(CTCF_Bone_marrow),GRanges(CTCF_Liver_epithelium))
hits0 <- findOverlaps(GRanges(CTCF_Bone_marrow[hits1@from,]),GRanges(CTCF_Liver_epithelium))


length(hits1@from)-length(hits0@from)
length(hits2@from)-length(hits0@from)
length(hits3@from)-length(hits0@from)

nrow(CTCF_Bone_marrow) - (length(hits1@from) + length(hits3@from) - length(hits0@from))
nrow(CTCF_B_Blood) - (length(hits1@from) + length(hits2@from) - length(hits0@from))
nrow(CTCF_Liver_epithelium) - (length(hits2@from) + length(hits3@from) - length(hits0@from))





hits1 <- findOverlaps(GRanges(Breast_epithelium_CTCF),GRanges(B_Blood_CTCF))

hits2 <- findOverlaps(GRanges(B_Blood_TF),GRanges(B_Blood_CTCF[hits1@to,]))

findOverlaps(GRanges(B_Blood_TF),GRanges(Breast_epithelium_CTCF))

library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
library(motifmatchr)
B_Blood_CTCF <- read.table('/media/shaoqizhu/easystore/CTCF/GSM935611_hg19_CTCF_B_Blood.narrowPeak')
B_Blood_Bcl11a <- read.table('/media/shaoqizhu/easystore/CTCF/GSM803388_hg19_Bcl11a_B_Blood.broadPeak')
B_Blood_Pax5 <- read.table('/media/shaoqizhu/easystore/CTCF/GSM803391_hg19_Pax5_B_Blood.broadPeak')
Breast_CTCF <- read.table('/media/shaoqizhu/easystore/CTCF/GSM822305_hg19_CTCF_Breast.narrowPeak')
Breast_FoxA1 <- read.table('/media/shaoqizhu/easystore/CTCF/GSM2257824_hg19_FOXA1_E3h_MCF7_peaks.bed')
colnames(B_Blood_CTCF)[1:3] <- c('chr','start','end')
colnames(B_Blood_Bcl11a)[1:3] <- c('chr','start','end')
colnames(B_Blood_Pax5)[1:3] <- c('chr','start','end')
colnames(Breast_CTCF)[1:3] <- c('chr','start','end')
colnames(Breast_FoxA1)[1:3] <- c('chr','start','end')

ctcf_motif <- matchMotifs(human_pwms_v1, GRanges(B_Blood_CTCF), genome = BSgenome.Hsapiens.UCSC.hg19)
match_num <-apply(ctcf_motif@assays@data$motifMatches,2,sum)
as.character(ctcf_motif@colData$name)[order(match_num,decreasing = T)[1:20]]
ctcf_idx <- which(ctcf_motif@assays@data$motifMatches[,which(ctcf_motif@colData$name=='CTCF')])
B_Blood_CTCF[ctcf_idx,]

ctcf_Breast_motif <- matchMotifs(human_pwms_v1, GRanges(Breast_CTCF), genome = BSgenome.Hsapiens.UCSC.hg19)
match_num <-apply(ctcf_Breast_motif@assays@data$motifMatches,2,sum)
as.character(ctcf_Breast_motif@colData$name)[order(match_num,decreasing = T)[1:20]]
ctcf_Breast_idx <- which(ctcf_Breast_motif@assays@data$motifMatches[,which(ctcf_Breast_motif@colData$name=='CTCF')])
Breast_CTCF[ctcf_Breast_idx,]

findOverlaps(GRanges(B_Blood_Bcl11a),GRanges(B_Blood_CTCF))
findOverlaps(GRanges(B_Blood_Bcl11a),GRanges(B_Blood_CTCF[-ctcf_idx,]))

findOverlaps(GRanges(B_Blood_Pax5),GRanges(B_Blood_CTCF))
findOverlaps(GRanges(B_Blood_Pax5),GRanges(B_Blood_CTCF[ctcf_idx,]))

hits1 <- findOverlaps(GRanges(Breast_CTCF[-ctcf_Breast_idx,]),GRanges(B_Blood_CTCF[-ctcf_idx,]))
hits2 <- findOverlaps(GRanges(Breast_CTCF[ctcf_Breast_idx,]),GRanges(B_Blood_CTCF[ctcf_idx,]))

findOverlaps(GRanges(B_Blood_Pax5),GRanges(B_Blood_CTCF[ctcf_idx,][hits2@to,]))
findOverlaps(GRanges(B_Blood_Pax5),GRanges(Breast_CTCF[ctcf_Breast_idx,][-hits2@from,]))

findOverlaps(GRanges(B_Blood_Pax5),GRanges(B_Blood_CTCF[-ctcf_idx,][hits1@to,]))
findOverlaps(GRanges(B_Blood_Pax5),GRanges(Breast_CTCF[-ctcf_Breast_idx,][-hits1@from,]))
findOverlaps(GRanges(B_Blood_Pax5),GRanges(Breast_CTCF))

findOverlaps(GRanges(Breast_FoxA1),GRanges(Breast_CTCF[-ctcf_Breast_idx,][hits1@from,]))
findOverlaps(GRanges(Breast_FoxA1),GRanges(B_Blood_CTCF[-ctcf_idx,][-hits1@to,]))


