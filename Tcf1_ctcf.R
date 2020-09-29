ctcf <- read.table('/media/shaoqizhu/easystore/CD8-HP/CTCF/macs2/CTCF_TCell_peaks.narrowPeak')
ctcf_summit <- read.table('/media/shaoqizhu/easystore/CD8-HP/CTCF/macs2/CTCF_TCell_summits.bed')
tcf1 <- read.table('/media/shaoqizhu/easystore/CD8-HP/CD8_HP_naive/NaiveCD8_TCF1_BAMPE_peaks_filtered.bed')
tcf1_summit <- read.table('/media/shaoqizhu/easystore/CD8-HP/CD8_HP_naive/19042_Feature_Tcf1_WT_CD8_peaks_summit.bed')
names(ctcf)[1:3] <- c('chr','start','end')
names(tcf1)[1:3] <- c('chr','start','end')
ctcf <- ctcf[ctcf$V9>-log10(0.05)&ctcf$V7>4,]
ctcf_summit <- ctcf_summit[rownames(ctcf),]

blacklist <- read.csv('/media/shaoqizhu/easystore/blacklist/mm9-blacklist.bed',sep = '\t',header = F)
names(blacklist)[1:3] <- c('chr','start','end')
hits_blacklist <- findOverlaps(GRanges(tcf1),GRanges(blacklist))
tcf1 <- tcf1[-hits_blacklist@from,]
tcf1_summit <- tcf1_summit[-hits_blacklist@from,]
hits_blacklist <- findOverlaps(GRanges(ctcf),GRanges(blacklist))
ctcf <- ctcf[-hits_blacklist@from,]
ctcf_summit <- ctcf_summit[-hits_blacklist@from,]

tcf1[order(tcf1$V9,decreasing = T)[1:3000],]
write.table(tcf1_summit[order(tcf1$V9,decreasing = T)[1:3000],],'/media/shaoqizhu/easystore/CD8-HP/homer/tcf1_qval_top3000.bed',sep = '\t',quote = F,col.names = F,row.names = F)

hits <- findOverlaps(GRanges(diff_peak_kmeans),GRanges(ctcf))
hits2 <- findOverlaps(GRanges(diff_peak_kmeans[hits@from,]),GRanges(tcf1))
cbind(1:7,table(diff_peak_kmeans[hits@from,][hits2@from,]$cluster_idx))

cbind(table(diff_peak_kmeans[hits@from,]$cluster_idx),
      table(diff_peak_kmeans$cluster_idx),
      table(diff_peak_kmeans[hits@from,]$cluster_idx)/table(diff_peak_kmeans$cluster_idx))
write.csv(diff_peak_kmeans[hits@from,],'~/Dropbox/CTCF_DNase.csv')

# tcf1 and ctcf chipseq overlap
hits <- findOverlaps(GRanges(ctcf),GRanges(tcf1))
head(tcf1[hits@to,])
hits
nrow(ctcf);nrow(tcf1);nrow(ctcf[hits@from,]); nrow(ctcf[-hits@from,]);nrow(tcf1[-hits@to,])

write.table(tcf1_summit[hits@to,][order(tcf1[hits@to,]$V9,decreasing = T),],
            '/media/shaoqizhu/easystore/CD8-HP/homer/region/tcf1_and_ctcf_ranked.bed',sep = '\t',quote = F,col.names = F,row.names = F)
write.table(tcf1_summit[hits@to,][order(tcf1[hits@to,]$V9,decreasing = T)[1:3000],],
            '/media/shaoqizhu/easystore/CD8-HP/homer/tcf1_and_ctcf_top3000.bed',sep = '\t',quote = F,col.names = F,row.names = F)

write.table(tcf1_summit[hits@to,],'/media/shaoqizhu/easystore/CD8-HP/homer/tcf1_and_ctcf.bed',sep = '\t',quote = F,col.names = F,row.names = F)
write.table(tcf1_summit[-hits@to,],'/media/shaoqizhu/easystore/CD8-HP/homer/tcf1_no_ctcf.bed',sep = '\t',quote = F,col.names = F,row.names = F)
write.table(ctcf_summit[hits@from,],'/media/shaoqizhu/easystore/CD8-HP/homer/ctcf_and_tcf1.bed',sep = '\t',quote = F,col.names = F,row.names = F)
write.table(ctcf_summit[-hits@from,],'/media/shaoqizhu/easystore/CD8-HP/homer/ctcf_no_tcf1.bed',sep = '\t',quote = F,col.names = F,row.names = F)

nrow(ctcf[-hits@from,])
tcf1[-hits@to,]

# tcf1 and ctcf chipseq overlap with promoters
hits <- findOverlaps(GRanges(ctcf),GRanges(tcf1))
hits1 <- findOverlaps(GRanges(tcf1[hits@to,]),GRanges(mm9_promoter))
hits2 <- findOverlaps(GRanges(ctcf[hits@from,]),GRanges(mm9_promoter))
nrow(tcf1[hits@to,][hits1@from,])
nrow(ctcf[hits@from,][hits2@from,])

tcf1_ctcf_prom_gene <- cbind(tcf1[hits@to,][hits1@from,],gene=mm9_promoter[hits1@to,]$gene_name)
tcf1_ctcf_prom_gene <- tcf1_ctcf_prom_gene[order(tcf1_ctcf_prom_gene$gene),]

RNA_seq_0 <- as.data.frame(RNA_seq[RNA_seq$gene_id%in%tcf1_ctcf_prom_gene$gene,])
row.names(RNA_seq_0) <- RNA_seq_0$gene_id
rna <- RNA_seq_0[as.character(tcf1_ctcf_prom_gene$gene),]
plot(log(tcf1_ctcf_prom_gene$V9),log(rna[,2]),xlab='tcf1 peak log(-log10(qval))',main='naive WT log2(RNA)')
plot(log(tcf1_ctcf_prom_gene$V9),log(rna[,5]/rna[,2]),xlab='tcf1 peak log(-log10(qval))',main='naive dKO/WT log2(RNA)')
examples <- which(log(rna[,5]/rna[,2])<(-1) & log(rna[,5]/rna[,2])>-10)
tcf1_ctcf_prom_gene[2105,]
rna[2105,1:10]
tcf1_ctcf_prom_gene[examples,]
rna[examples,1:10]

WT_na_dKO_na <- read.csv('/media/shaoqizhu/easystore/CD8-HP/diff_peaks/WT_na-dKO_na.bed',sep = '\t',header = F)
names(WT_na_dKO_na)[1:3] <- c('chr','start','end')
hits3 <- findOverlaps(GRanges(tcf1[hits@to,][-hits1@from,]),GRanges(WT_na_dKO_na[WT_na_dKO_na$V7<0,]))
write.table(tcf1_summit[hits@to,][-hits1@from,][hits3@from,],'/media/shaoqizhu/easystore/CD8-HP/homer/region/tcf1_and_ctcf_dist_na_WT_mt_dKO.bed',sep = '\t',quote = F,col.names = F,row.names = F)
hits4 <- findOverlaps(GRanges(tcf1[hits@to,][-hits1@from,]),GRanges(WT_na_dKO_na[WT_na_dKO_na$V7>0,]))
write.table(tcf1_summit[hits@to,][-hits1@from,][hits4@from,],'/media/shaoqizhu/easystore/CD8-HP/homer/region/tcf1_and_ctcf_dist_na_WT_lt_dKO.bed',sep = '\t',quote = F,col.names = F,row.names = F)
hits5 <- findOverlaps(GRanges(tcf1[hits@to,][-hits1@from,][hits3@from,]),GRanges(mm9_deg_50kb))
hits6 <- findOverlaps(GRanges(tcf1[hits@to,][-hits1@from,][hits4@from,]),GRanges(mm9_deg_50kb))
table(mm9_deg_50kb[hits6@to,]$cluster)



write.table(tcf1[hits@to,][hits1@from,],'/media/shaoqizhu/easystore/CD8-HP/homer/region/tcf1_and_ctcf_prom_peak.bed',sep = '\t',quote = F,col.names = F,row.names = F)
tcf1_and_ctcf_prom_na_dKO <- read.csv('/media/shaoqizhu/easystore/CD8-HP/DNase/tcf1_and_ctcf_prom_na_dKO.bed',sep = '\t',header = F)
tcf1_and_ctcf_prom_na_WT <- read.csv('/media/shaoqizhu/easystore/CD8-HP/DNase/tcf1_and_ctcf_prom_na_WT.bed',sep = '\t',header = F)
plot(log2(tcf1_and_ctcf_prom_na_WT$V9),log2(tcf1_and_ctcf_prom_na_WT$V11),xlab='tcf1 peak -log10(qval)',ylab='naive WT log2(coverage)')
plot(log2(tcf1_and_ctcf_prom_na_WT$V9),log2(1831553/1173094*tcf1_and_ctcf_prom_na_dKO$V11/tcf1_and_ctcf_prom_na_WT$V11),xlab='tcf1 peak -log10(qval)',ylab='naive dKO/WT log2(coverage)')

sum(tcf1_and_ctcf_prom_na_dKO$V11)
sum(tcf1_and_ctcf_prom_na_WT$V11)

write.table(tcf1_summit[hits@to,][hits1@from,],'/media/shaoqizhu/easystore/CD8-HP/homer/region/tcf1_and_ctcf_prom.bed',sep = '\t',quote = F,col.names = F,row.names = F)
write.table(tcf1_summit[hits@to,][-hits1@from,],'/media/shaoqizhu/easystore/CD8-HP/homer/region/tcf1_and_ctcf_dist.bed',sep = '\t',quote = F,col.names = F,row.names = F)
write.table(ctcf_summit[hits@from,][hits2@from,],'/media/shaoqizhu/easystore/CD8-HP/homer/region/ctcf_and_tcf1_prom.bed',sep = '\t',quote = F,col.names = F,row.names = F)
write.table(ctcf_summit[hits@from,][-hits2@from,],'/media/shaoqizhu/easystore/CD8-HP/homer/region/ctcf_and_tcf1_dist.bed',sep = '\t',quote = F,col.names = F,row.names = F)

write.table(tcf1_summit[hits@to,][hits1@from,][order(tcf1[hits@to,][hits1@from,]$V9,decreasing = T),],'/media/shaoqizhu/easystore/CD8-HP/homer/region/tcf1_and_ctcf_prom_ranked.bed',sep = '\t',quote = F,col.names = F,row.names = F)
write.table(tcf1_summit[hits@to,][-hits1@from,][order(tcf1[hits@to,][-hits1@from,]$V9,decreasing = T),],'/media/shaoqizhu/easystore/CD8-HP/homer/region/tcf1_and_ctcf_dist_ranked.bed',sep = '\t',quote = F,col.names = F,row.names = F)
head(tcf1[hits@to,][-hits1@from,][order(tcf1[hits@to,][-hits1@from,]$V9,decreasing = T),])

write.table(tcf1_summit[hits@to,][hits1@from,][order(tcf1[hits@to,][hits1@from,]$V9,decreasing = T)[1:2000],],
            '/media/shaoqizhu/easystore/CD8-HP/homer/region/tcf1_and_ctcf_prom_top2000.bed',sep = '\t',quote = F,col.names = F,row.names = F)
write.table(ctcf_summit[hits@from,][hits2@from,][order(ctcf[hits@from,][hits2@from,]$V9,decreasing = T)[1:2000],],
            '/media/shaoqizhu/easystore/CD8-HP/homer/region/ctcf_and_tcf1_prom_top2000.bed',sep = '\t',quote = F,col.names = F,row.names = F)
tcf1_ctcf_2000 <- tcf1[hits@to,][hits1@from,][order(tcf1[hits@to,][hits1@from,]$V9,decreasing = T)[1:2000],]
ctcf_tcf1_2000 <- ctcf[hits@from,][hits2@from,][order(ctcf[hits@from,][hits2@from,]$V9,decreasing = T)[1:2000],]
hits3 <- findOverlaps(GRanges(tcf1_ctcf_2000),GRanges(ctcf_tcf1_2000))
write.table(tcf1_summit[hits@to,][hits1@from,][order(tcf1[hits@to,][hits1@from,]$V9,decreasing = T)[1:2000],][hits3@from,],
            '/media/shaoqizhu/easystore/CD8-HP/homer/region/tcf1_and_ctcf_prom_top2000_overlap1083.bed',sep = '\t',quote = F,col.names = F,row.names = F)

hits3 <- findOverlaps(GRanges(tcf1[hits@to,]),GRanges(stat5_wt_24h))
#hits4 <- findOverlaps(GRanges(tcf1[hits@to,][hits3@from,]),GRanges(stat5_dko_24h))
hits0 <- findOverlaps(GRanges(tcf1[hits@to,][hits3@from,]),GRanges(diff_peak_kmeans))
diff_peak_kmeans[hits0@to,][201:400,]


hits1_plus <- findOverlaps(GRanges(tcf1[hits@to,]),GRanges(mm9_promoter[mm9_promoter$strand=='+',]))
hits1_minus <- findOverlaps(GRanges(tcf1[hits@to,]),GRanges(mm9_promoter[mm9_promoter$strand=='-',]))
hits2_plus <- findOverlaps(GRanges(ctcf[hits@from,]),GRanges(mm9_promoter[mm9_promoter$strand=='+',]))
hits2_minus <- findOverlaps(GRanges(ctcf[hits@from,]),GRanges(mm9_promoter[mm9_promoter$strand=='-',]))

hits_plus <- findOverlaps(GRanges(ctcf),GRanges(mm9_promoter[mm9_promoter$strand=='+',]))
hits_minus <- findOverlaps(GRanges(ctcf),GRanges(mm9_promoter[mm9_promoter$strand=='-',]))
nrow(ctcf[hits_plus@from,])
nrow(ctcf[hits_minus@from,])

tcf1[hits@to,][hits1_minus@from,]
ctcf[hits@from,][hits2_minus@from,]
write.table(tcf1_summit[hits@to,][hits1_plus@from,],'/media/shaoqizhu/easystore/CD8-HP/homer/region/tcf1_and_ctcf_prom_plus.bed',sep = '\t',quote = F,col.names = F,row.names = F)
write.table(tcf1_summit[hits@to,][hits1_minus@from,],'/media/shaoqizhu/easystore/CD8-HP/homer/region/tcf1_and_ctcf_prom_minus.bed',sep = '\t',quote = F,col.names = F,row.names = F)
write.table(ctcf_summit[hits@from,][hits2_plus@from,],'/media/shaoqizhu/easystore/CD8-HP/homer/region/ctcf_and_tcf1_prom_plus.bed',sep = '\t',quote = F,col.names = F,row.names = F)
write.table(ctcf_summit[hits@from,][hits2_minus@from,],'/media/shaoqizhu/easystore/CD8-HP/homer/region/ctcf_and_tcf1_prom_minus.bed',sep = '\t',quote = F,col.names = F,row.names = F)

hits_plus <- findOverlaps(GRanges(tcf1),GRanges(mm9_promoter[mm9_promoter$strand=='+',]))
hits_minus <- findOverlaps(GRanges(tcf1),GRanges(mm9_promoter[mm9_promoter$strand=='-',]))
length(hits_plus@from)
length(hits_minus@from)

tcf1_and_ctcf.bed
motif_output <- read.csv('/media/shaoqizhu/easystore/CD8-HP/homer/MotifOutput/outputfile.txt',sep = '\t')
motif_output <- motif_output[order(motif_output$PositionID),]
motif_copeak <- motif_output[motif_output$PositionID %in% names(which(table(motif_output$PositionID)>=2)),]
table(motif_output$Motif.Name)
agg <- aggregate(motif_copeak$Motif.Name,list(motif_copeak$PositionID),function(x){length(unique(x))})
motif_cobind <-motif_copeak[motif_copeak$PositionID %in% agg$Group.1[agg$x==2],]
motif_cobind <- motif_cobind[order(motif_cobind$PositionID,motif_cobind$Motif.Name),]
hist(aggregate(motif_cobind$Offset,list(motif_cobind$PositionID),function(x) x[2]-x[1])$x,
     xlab = 'Tcf1-Ctcf',main='Tcf1-Ctcf motif distance')

table(aggregate(motif_cobind$Strand,list(motif_cobind$PositionID),function(x){length(unique(x))})$x)

names(which(table(motif_output$PositionID)>=2))

hist(motif_output$Offset[motif_output$Motif.Name=='Tcf7'])

ctcf_tcf1 <- as.data.frame(reduce(GRanges(rbind(ctcf[hits@from,1:3],tcf1[hits@to,1:3]))))[1:3]

hits <- findOverlaps(GRanges(ctcf_tcf1),GRanges(mm9_deg_50kb))

table(mm9_deg_50kb[hits@to,"cluster"])
cbind(table(mm9_deg$cluster),
      table(aggregate(mm9_deg_50kb[hits@to,"cluster"],list(mm9_deg_50kb[hits@to,"gene_name"]),unique)$x))

write.table(ctcf_tcf1,'/media/shaoqizhu/easystore/CD8-HP/ctcf_tcf1.bed',sep = '\t',quote = F,col.names = F,row.names = F)
write.table(diff_peak_kmeans[diff_peak_kmeans$cluster_idx==1,],'/media/shaoqizhu/easystore/CD8-HP/diff_peak_1.bed',sep = '\t',quote = F,col.names = F,row.names = F)

ctcf_tcf1_spamo <- read.table('/media/shaoqizhu/easystore/CD8-HP/meme/ctcf_tcf1_spamo.bed')
names(ctcf_tcf1_spamo) <- c('chr','start','end')
hits <- findOverlaps(GRanges(diff_peak_kmeans),GRanges(ctcf_tcf1_spamo))
diff_peak_kmeans[hits@from,]

motif_ix <- matchMotifs(mouse_pwms_v1, GRanges(tcf1[order(tcf1$V5,decreasing = T)[1:1000],1:3]),
                        genome = BSgenome.Mmusculus.UCSC.mm9)

match_num <- apply(motif_ix@assays@data$motifMatches,2,sum)

motif_match <- motif_ix@colData[order(match_num,decreasing = T)[1:100],]

as.data.frame(cbind(as.character(motif_match),sort(match_num,decreasing = T)[1:100]))

which(as.character(motif_match)=='Lef')
grep('Tcf',as.character(motif_match))


CD8_Runx3 <- read.table('/media/shaoqizhu/easystore/CD8-HP/CTCF/TF_CD8/GSE50130_Resting_CD8_Runx3_peaks_to_small.bed')
CD8_CBF <- read.table('/media/shaoqizhu/easystore/CD8-HP/CTCF/TF_CD8/GSM2471915_Naive_CD8_CBF_peaks.bed')
colnames(CD8_Runx3)[1:3] <- c('chr','start','end')
colnames(CD8_CBF)[1:3] <- c('chr','start','end')

findOverlaps(GRanges(tcf1),GRanges(CD8_Runx3))
findOverlaps(GRanges(tcf1),GRanges(CD8_CBF))
findOverlaps(GRanges(tcf1),GRanges(stat5_wt_24h))
hits <- findOverlaps(GRanges(CD8_CBF),GRanges(CD8_Runx3))
findOverlaps(GRanges(ctcf),GRanges(CD8_CBF[hits@from,]))

library(motifmatchr)
ctcf <- ctcf[-which(ctcf$chr=='chr1_random'),]
ctcf_motif <- matchMotifs(mouse_pwms_v1, GRanges(ctcf), genome = BSgenome.Mmusculus.UCSC.mm9)
match_num <-apply(ctcf_motif@assays@data$motifMatches,2,sum)
as.character(ctcf_motif@colData$name)[order(match_num,decreasing = T)[1:20]]
ctcf_idx <- which(ctcf_motif@assays@data$motifMatches[,which(ctcf_motif@colData$name=='Ctcf')])

findOverlaps(GRanges(stat5_wt_24h),GRanges(ctcf[ctcf_idx,]))

hits <- findOverlaps(GRanges(tcf1),GRanges(ctcf))
tcf1_order <- tcf1[hits@from,]#[order(tcf1[hits@from,]$V6,decreasing = T)[1:500],]
hits1 <- findOverlaps(GRanges(mm9_promoter),GRanges(tcf1_order))

hits2 <- findOverlaps(GRanges(tcf1_order[hits1@to,]),GRanges(merge))
ratio_prom <- log(RPKM[hits2@to,]$`12_dKO_na`/RPKM[hits2@to,]$`11_WT_na`)
hist(ratio_prom[ratio_prom>-3 & ratio_prom<3],breaks = seq(-3,3,0.3),main='Naive log(dko/wt) DNase on promoter')

hits2 <- findOverlaps(GRanges(tcf1_order[-hits1@to,]),GRanges(merge))
ratio_enhancer <- log(RPKM[hits2@to,]$`12_dKO_na`/RPKM[hits2@to,]$`11_WT_na`)
hist(ratio_enhancer[ratio_enhancer>-3 & ratio_enhancer<3],breaks = seq(-3,3,0.3),main='Naive dko/wt DNase on enhancer')

df <- data.frame(
  x = c(ratio_prom[ratio_prom>-3 & ratio_prom<3],
        ratio_enhancer[ratio_enhancer>-3 & ratio_enhancer<3]),
  Tcf1_coverage = factor(c(rep('promoter',length(ratio_prom[ratio_prom>-3 & ratio_prom<3])),
                           rep('enhancer',length(ratio_enhancer[ratio_enhancer>-3 & ratio_enhancer<3]))))
)
ggplot(df, aes(x,colour=Tcf1_coverage)) + stat_ecdf(geom = "step") + xlab('log(DNase DKO/WT)')

hits <- findOverlaps(GRanges(mm9_promoter),GRanges(tcf1))
tcf1_order <- tcf1[hits@to,]
hits1 <- findOverlaps(GRanges(tcf1_order),GRanges(ctcf))

hits2 <- findOverlaps(GRanges(tcf1_order[hits1@from,]),GRanges(merge))
ratio_ctcf <- log(RPKM[hits2@to,]$`12_dKO_na`/RPKM[hits2@to,]$`11_WT_na`)
hist(ratio_ctcf[ratio_ctcf>-3 & ratio_ctcf<3],breaks = seq(-3,3,0.3),main='Naive log(dko/wt) DNase on Tcf1+CTCF+ promoter')

hits2 <- findOverlaps(GRanges(tcf1_order[-hits1@from,]),GRanges(merge))
ratio_noctcf <- log(RPKM[hits2@to,]$`12_dKO_na`/RPKM[hits2@to,]$`11_WT_na`)
hist(ratio_noctcf[ratio_noctcf>-3 & ratio_noctcf<3],breaks = seq(-3,3,0.3),main='Naive log(dko/wt) DNase on Tcf1+CTCF- promoter')

df <- data.frame(
  x = c(ratio_ctcf[ratio_ctcf>-3 & ratio_ctcf<3],
        ratio_noctcf[ratio_noctcf>-3 & ratio_noctcf<3]),
  Tcf1_coverage = factor(c(rep('Tcf1+CTCF+',length(ratio_ctcf[ratio_ctcf>-3 & ratio_ctcf<3])),
                           rep('Tcf1+CTCF-',length(ratio_noctcf[ratio_noctcf>-3 & ratio_noctcf<3]))))
)
ggplot(df, aes(x,colour=Tcf1_coverage)) + stat_ecdf(geom = "step") + xlab('log(DNase DKO/WT)')




write.table(tcf1_order[-hits1@to,c(1,2,3,10)],'/media/shaoqizhu/easystore/tcf1_ctcf_enhancer.bed',sep = '\t',quote = F,col.names = F,row.names = F)
write.table(tcf1_order[hits1@to,c(1,2,3,10)],'/media/shaoqizhu/easystore/tcf1_ctcf_promoter.bed',sep = '\t',quote = F,col.names = F,row.names = F)
write.table(tcf1_order[-hits1@to,c(1,2,3,10)],'/media/shaoqizhu/easystore/tcf1_ctcf_motif_enhancer.bed',sep = '\t',quote = F,col.names = F,row.names = F)
write.table(tcf1_order[-hits1@to,c(1,2,3,10)],'/media/shaoqizhu/easystore/tcf1_noctcf_enhancer.bed',sep = '\t',quote = F,col.names = F,row.names = F)
write.table(tcf1[-hits1@to,c(1,2,3,10)],'/media/shaoqizhu/easystore/tcf1_enhancer.bed',sep = '\t',quote = F,col.names = F,row.names = F)

tcf1_ctcf_nomotif_enhancer <- read.csv('/media/shaoqizhu/easystore/CTCF/tcf1/tcf1_ctcf_nomotif_enhancer.txt',sep='\t')
tcf1_ctcf_motif_enhancer <- read.csv('/media/shaoqizhu/easystore/CTCF/tcf1/tcf1_ctcf_motif_enhancer.txt',sep='\t')
tcf1_noctcf_enhancer <- read.csv('/media/shaoqizhu/easystore/CTCF/tcf1/tcf1_noctcf_enhancer.txt',sep='\t')
tcf1_enhancer <- read.csv('/media/shaoqizhu/easystore/CTCF/tcf1/tcf1_enhancer.txt',sep='\t')

for (i in 1:8){
  print(length(intersect(tcf1_enhancer$X..GREAT.version.4.0.4,mm9_deg$gene_name[mm9_deg$cluster==i])))
}
table(mm9_deg$cluster)

for (i in 1:8){
  print(length(intersect(mm9_promoter[hits1@from,]$gene_name,mm9_deg$gene_name[mm9_deg$cluster==i])))
}

tcf1_gene = mm9_promoter[hits1@from,]
tcf1_gene$wt_rpkm=0
tcf1_gene$dko_rpkm=0

for (i in 1:nrow(tcf1_gene)){
  tcf1_gene$wt_rpkm[i] = RNA_seq$WT_0h_5n[which(RNA_seq$gene_id==tcf1_gene$gene_name[i])]
  tcf1_gene$dko_rpkm[i] = RNA_seq$DKO_0h_5n[which(RNA_seq$gene_id==tcf1_gene$gene_name[i])]
}

plot(log(tcf1_order$V9[hits1@to]),log2(tcf1_gene$dko_rpkm/tcf1_gene$wt_rpkm))

write.csv(tcf1_gene$gene_name[which(log2(tcf1_gene$dko_rpkm/tcf1_gene$wt_rpkm)< -1)],
          '/media/shaoqizhu/easystore/gene.csv')

hits <- findOverlaps(GRanges(ctcf),GRanges(tcf1))
findOverlaps(GRanges(ctcf[hits@from,]),GRanges(mm9_promoter))




