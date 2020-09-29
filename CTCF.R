library(rtracklayer)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(motifmatchr)
load("~/Downloads/chromVARmotifs-master/data/human_pwms_v1.rda")

gtf <- import('/media/shaoqizhu/easystore/Annotation/hg19_2015.gtf')
gtf <- as.data.frame(gtf)

gene.name = sort(unique(gtf$gene_name))
gtf_reduced <- as.data.frame(matrix(NA,ncol = 6,nrow = 0))
for (x in 1:length(gene.name)){
  gene = gtf[gtf$gene_name==gene.name[x],]
  transcript_chr = aggregate(gene$seqnames,list(gene$transcript_id),unique)[,2]
  transcript_start = aggregate(gene$start,list(gene$transcript_id),min)[,2]
  transcript_end = aggregate(gene$end,list(gene$transcript_id),max)[,2]
  transcript_strand = as.character(aggregate(gene$strand,list(gene$transcript_id),unique)[,2])
  transcript_ranges <- as.data.frame(reduce(GRanges(transcript_chr,IRanges(transcript_start,transcript_end),transcript_strand)))
  gtf_reduced = rbind(gtf_reduced,cbind(transcript_ranges,gene_name=rep(gene.name[x],nrow(transcript_ranges))))
}
hg19_reduced <- gtf_reduced
hg19_reduced <- hg19_reduced[c(1:4,6,5)]
write.table(hg19_reduced,'/media/shaoqizhu/easystore/Annotation/hg19_2015_reduced.bed',quote = F,col.names = F,sep='\t',row.names = F)
hg19_reduced <- read.table('/media/shaoqizhu/easystore/Annotation/hg19_2015_reduced.bed')
colnames(hg19_reduced) <- colnames(mm9_reduced)
hg19_promoter <- hg19_reduced
for (i in 1:nrow(hg19_reduced)){
  if(hg19_reduced$strand[i]=='+'){
    ifelse(hg19_reduced$start[i]>1000,
           hg19_promoter$start[i] <- hg19_reduced$start[i]-1000,
           hg19_promoter$start[i] <- 0)
    hg19_promoter$end[i] <- hg19_reduced$start[i]+1000
  }
  if(hg19_reduced$strand[i]=='-'){
    ifelse(hg19_reduced$end[i]>1000,
           hg19_promoter$start[i] <- hg19_reduced$end[i]-1000,
           hg19_promoter$start[i] <- 0)
    hg19_promoter$end[i] <- hg19_reduced$end[i]+1000
  }
}

#============================ CTCF chipseq
Venn <- function(S1,S2,S3){
  hits1 <- findOverlaps(GRanges(S1),GRanges(S2))
  hits2 <- findOverlaps(GRanges(S1),GRanges(S3))
  hits3 <- findOverlaps(GRanges(S2),GRanges(S3))
  hits0 <- findOverlaps(GRanges(S1[unique(hits1@from),]),GRanges(S3))
  A=nrow(S1);B=nrow(S2);D=nrow(S3);
  c <- length(unique(hits1@from))-length(unique(hits0@from))
  e <- length(unique(hits2@from))-length(unique(hits0@from))
  f <- length(unique(hits3@from))-length(unique(hits0@from))
  g <- length(unique(hits0@from))
  a=A-c-e-g;b=B-c-f-g;d=D-e-f-g
  paste(a,b,c,d,e,f,g,A,B,D,sep = ',')
}

PAX5_GM12878 <- read.table('/media/shaoqizhu/easystore/CTCF/ChIA-PET/GM12878/PAX5_GM12878_ENCFF969EMZ.bed')
CEBPB_Hela <- read.table('/media/shaoqizhu/easystore/CTCF/ChIA-PET/HeLa/CEBPB_HeLa_ENCFF002CSA.bed')
GATA1_K562 <- read.table('/media/shaoqizhu/easystore/CTCF/ChIA-PET/K562/GATA1_K562_ENCFF178NBS.bed')
colnames(PAX5_GM12878)[1:3] <- c('chr','start','end')
colnames(CEBPB_Hela)[1:3] <- c('chr','start','end')
colnames(GATA1_K562)[1:3] <- c('chr','start','end')

CTCF_GM12878 <- read.table('/media/shaoqizhu/easystore/CTCF/ChIA-PET/CTCF_GM12878_ENCFF000ARG-rmdup-W200-G200.scoreisland')
CTCF_Hela <- read.table('/media/shaoqizhu/easystore/CTCF/ChIA-PET/CTCF_HeLa_ENCFF000BAJ_rmdup-W200-G200.scoreisland')
CTCF_K562 <- read.table('/media/shaoqizhu/easystore/CTCF/ChIA-PET/CTCF_K562_ENCFF000YLT_rmdup-W200-G200.scoreisland')
colnames(CTCF_GM12878)[1:3] <- c('chr','start','end')
colnames(CTCF_Hela)[1:3] <- c('chr','start','end')
colnames(CTCF_K562)[1:3] <- c('chr','start','end')

CTCF_GM12878_prom <- CTCF_GM12878[findOverlaps(GRanges(CTCF_GM12878),GRanges(hg19_promoter))@from,]
CTCF_Hela_prom <- CTCF_Hela[findOverlaps(GRanges(CTCF_Hela),GRanges(hg19_promoter))@from,]
CTCF_K562_prom <- CTCF_K562[findOverlaps(GRanges(CTCF_K562),GRanges(hg19_promoter))@from,]
CTCF_GM12878_enhancer <- CTCF_GM12878[-findOverlaps(GRanges(CTCF_GM12878),GRanges(hg19_promoter))@from,]
CTCF_Hela_enhancer <- CTCF_Hela[-findOverlaps(GRanges(CTCF_Hela),GRanges(hg19_promoter))@from,]
CTCF_K562_enhancer <- CTCF_K562[-findOverlaps(GRanges(CTCF_K562),GRanges(hg19_promoter))@from,]

Venn(CTCF_GM12878_prom,CTCF_Hela_prom,CTCF_K562_prom)
Venn(CTCF_GM12878_enhancer,CTCF_Hela_enhancer,CTCF_K562_enhancer)

hits1 <- findOverlaps(GRanges(CTCF_Hela[-motif_idx_H,]),GRanges(CTCF_GM12878[-motif_idx_G,]))
hits2 <- findOverlaps(GRanges(CTCF_Hela[-motif_idx_H,][-unique(hits1@from),]),GRanges(CTCF_K562[-motif_idx_K,]))
CTCF_Hela_unique <- CTCF_Hela[-motif_idx_H,][-unique(hits1@from),][-unique(hits2@from),]
findOverlaps(GRanges(CEBPB_Hela),GRanges(CTCF_Hela_unique))

hits1 <- findOverlaps(GRanges(CTCF_GM12878[-motif_idx_G,]),GRanges(CTCF_Hela[-motif_idx_H,]))
hits2 <- findOverlaps(GRanges(CTCF_GM12878[-motif_idx_G,][-unique(hits1@from),]),GRanges(CTCF_K562[-motif_idx_K,]))
CTCF_GM12878_unique <- CTCF_GM12878[-motif_idx_G,][-unique(hits1@from),][-unique(hits2@from),]
hits0 <- findOverlaps(GRanges(BCL11A_GM12878),GRanges(CTCF_GM12878_unique))
findOverlaps(GRanges(PAX5_GM12878),GRanges(CTCF_GM12878_unique[hits0@to,]))

hits1 <- findOverlaps(GRanges(CTCF_K562[-motif_idx_K,]),GRanges(CTCF_Hela[-motif_idx_H,]))
hits2 <- findOverlaps(GRanges(CTCF_K562[-motif_idx_K,][-unique(hits1@from),]),GRanges(CTCF_GM12878[-motif_idx_G,]))
CTCF_K562_unique <- CTCF_K562[-motif_idx_K,][-unique(hits1@from),][-unique(hits2@from),]
findOverlaps(GRanges(GATA1_K562),GRanges(CTCF_K562_unique))

write.table(CTCF_K562_unique[order(CTCF_K562_unique$V4,decreasing = T)[1:5000],],
            '/media/shaoqizhu/easystore/CTCF/CTCF_K562_unique_top5000.bed',quote = F,col.names = F,row.names = F,sep = '\t')


write.table(CTCF_GM12878_unique,'/media/shaoqizhu/easystore/CTCF/CTCF_GM12878_unique.bed',quote = F,col.names = F,row.names = F,sep = '\t')
write.table(CTCF_Hela_unique,'/media/shaoqizhu/easystore/CTCF/CTCF_Hela_unique.bed',quote = F,col.names = F,row.names = F,sep = '\t')
write.table(CTCF_K562_unique,'/media/shaoqizhu/easystore/CTCF/CTCF_K562_unique.bed',quote = F,col.names = F,row.names = F,sep = '\t')


#================================== CHIA-PET

CHIA_GM12878 <- read.table('/media/shaoqizhu/easystore/CTCF/ChIA-PET/GM12878/GM12878_MICC_loops_IAB2FDR0.05.txt',header = T)
CHIA_Hela <- read.table('/media/shaoqizhu/easystore/CTCF/ChIA-PET/HeLa/HeLa_MICC_loops_IAB2FDR0.05.txt',header = T)
CHIA_K562 <- read.table('/media/shaoqizhu/easystore/CTCF/ChIA-PET/K562/K562_MICC_loops_IAB2FDR0.05.txt',header = T)
colnames(CHIA_GM12878)[1:3] <- c('chr','start','end');colnames(CHIA_GM12878)[4:6] <- c('chr','start','end')
colnames(CHIA_Hela)[1:3] <- c('chr','start','end');colnames(CHIA_Hela)[4:6] <- c('chr','start','end')
colnames(CHIA_K562)[1:3] <- c('chr','start','end');colnames(CHIA_K562)[4:6] <- c('chr','start','end')
CHIA_GM12878_1 <- CHIA_GM12878[,1:3];CHIA_GM12878_2 <- CHIA_GM12878[,4:6]
CHIA_Hela_1 <- CHIA_Hela[,1:3];CHIA_Hela_2 <- CHIA_Hela[,4:6]
CHIA_K562_1 <- CHIA_K562[,1:3];CHIA_K562_2 <- CHIA_K562[,4:6]

overlap_CHIA <- function(A,B){
  A1 <- A[,1:3]; A2 <- A[,4:6];B1 <- B[,1:3]; B2 <- B[,4:6]
  hits1 <- findOverlaps(GRanges(A1),GRanges(B1))
  hits2 <- findOverlaps(GRanges(A2),GRanges(B2))
  left <- paste(hits1@from,hits1@to,sep = ':')
  right <- paste(hits2@from,hits2@to,sep = ':')
  as.data.frame(cbind(hits1@from[which(left%in%right)],hits2@to[which(right%in%left)]))
}
library(GenomicRanges)
idx1 <- overlap_CHIA(CHIA_GM12878,CHIA_Hela)
idx2 <- overlap_CHIA(CHIA_GM12878,CHIA_K562)
idx3 <- overlap_CHIA(CHIA_K562,CHIA_Hela)

CHIA_GM12878[idx1$V1[200],]
CHIA_Hela[idx1$V2[200],]

idx4 <- overlap_CHIA(CHIA_GM12878[-unique(idx1$V1),],CHIA_K562)
nrow(CHIA_GM12878[-unique(idx1$V1),][-unique(idx4$V1),])

idx5 <- overlap_CHIA(CHIA_Hela[-unique(idx1$V2),],CHIA_K562)
nrow(CHIA_Hela[-unique(idx1$V2),][-unique(idx5$V1),])

idx6 <- overlap_CHIA(CHIA_K562[-unique(idx3$V1),],CHIA_GM12878)
CHIA_K562[-unique(idx3$V1),][-unique(idx6$V1),]


#===================================== Graph

reduce(GRanges(unique(rbind(CHIA_GM12878[,1:3],CHIA_GM12878[,1:3]))))
CHIA_GM_anchor <- rbind(CHIA_GM12878[,1:3],CHIA_GM12878[,4:6])
anchor_order <- CHIA_GM_anchor[order(CHIA_GM_anchor$chr,CHIA_GM_anchor$start,CHIA_GM_anchor$end),]

anchor_unique <- cbind(aggregate(anchor_order[,1:3],list(anchor_order$chr,anchor_order$start),unique)[,3:5],
                       num=aggregate(anchor_order$chr,list(anchor_order$chr,anchor_order$start),length)$x)      
anchor_unique <- anchor_unique[order(anchor_unique$chr,anchor_unique$start,anchor_unique$end),]
rownames(anchor_unique) <- 1:nrow(anchor_unique)

anchor_idx <- as.numeric()
for (i in 1:nrow(anchor_unique)){anchor_idx <- c(anchor_idx,rep(i,anchor_unique$num[i]))}
anchor_order$idx <- anchor_idx
anchor_edge <- as.data.frame(cbind(A = anchor_order[order(as.numeric(rownames(anchor_order))),][1:(nrow(anchor_order)/2),],
                     B = anchor_order[order(as.numeric(rownames(anchor_order))),][(nrow(anchor_order)/2+1):nrow(anchor_order),],
                     IAB = CHIA_GM12878$IAB))


CHIA_GM12878[which(anchor_edge$A.idx %in% g1_idx & anchor_edge$B.idx %in% g1_idx),]
anchor_edge[which(anchor_edge$A.idx %in% g1_idx & anchor_edge$B.idx %in% g1_idx),]

hits <- findOverlaps(GRanges(PAX5_GM12878),GRanges(CTCF_GM12878))
hits1 <- findOverlaps(GRanges(PAX5_GM12878),GRanges(anchor_unique))
unique(hits1@to)

nrow(CHIA_GM12878[which(anchor_edge$A.idx %in% unique(hits1@to)
                   & anchor_edge$B.idx %in% unique(hits1@to)),])

random_anchor <- sample(1:nrow(anchor_unique),4163)
nrow(CHIA_GM12878[which(anchor_edge$A.idx %in% random_anchor
                        & anchor_edge$B.idx %in% random_anchor),])


anchor_cluster <- as.data.frame(cbind(as.numeric(g_cluster$names),g_cluster$membership))
PAX5_anchor <- sort(table(anchor_cluster$V2[anchor_cluster$V1 %in% unique(hits@to)]),decreasing = T)[1:10]

table(g_cluster$membership)[names(table(g_cluster$membership)) %in% names(PAX5_anchor)]
mm9_deg_pro



library(igraph)
cluster_louvain(g)

g <- graph.data.frame(anchor_edge[,c(4,8,9)],directed = F)

g_cluster <- cluster_louvain(g)
g1_idx = as.numeric(g_cluster$names[g_cluster$membership==8])

which(anchor_edge$A.idx %in% g1_idx & anchor_edge$B.idx %in% g1_idx)
g1 <- graph.data.frame(anchor_edge[anchor_edge$A.idx %in% g1_idx & anchor_edge$B.idx %in% g1_idx,c(4,8,9)],
                       directed = T)

plot(g1,layout=layout_with_gem, vertex.size=15, vertex.shapes='nil',vertex.color="green")


#===================================== CHIA-PET CTCF
hits_G1 <- findOverlaps(GRanges(CTCF_GM12878),GRanges(CHIA_GM12878_1))
hits_G2 <- findOverlaps(GRanges(CTCF_GM12878),GRanges(CHIA_GM12878_2))
unique(hits_G2@to)
loop_with_ctcf <- intersect(unique(hits_G1@to),unique(hits_G2@to))

CTCF_GM12878[hits_G1@from[which(hits_G1@to==loop_with_ctcf[2])],]
CTCF_GM12878[hits_G2@from[which(hits_G2@to==loop_with_ctcf[2])],]
boxplot(log(CHIA_GM12878[loop_with_ctcf,]$IAB),log(CHIA_GM12878[-loop_with_ctcf,]$IAB))

#=================================TF and CHIA-PET


hits1 <- findOverlaps(GRanges(PAX5_GM12878),GRanges(CTCF_GM12878_prom))
hits2 <- findOverlaps(GRanges(PAX5_GM12878[hits1@from,]),GRanges(CHIA_GM12878_1))
hits3 <- findOverlaps(GRanges(PAX5_GM12878),GRanges(CTCF_GM12878_enhancer))
hits4 <- findOverlaps(GRanges(PAX5_GM12878[hits3@from,]),GRanges(CHIA_GM12878_1))
boxplot(log(CHIA_GM12878$IAB[unique(hits2@to)]),log(CHIA_GM12878$IAB[unique(hits4@to)]),names=c('promoter','enhancer'))
ks.test(log(CHIA_GM12878$IAB[unique(hits2@to)]),log(CHIA_GM12878$IAB[unique(hits4@to)]))

hits1 <- findOverlaps(GRanges(PAX5_GM12878),GRanges(CTCF_GM12878))
hits2 <- findOverlaps(GRanges(PAX5_GM12878[hits1@from,]),GRanges(CHIA_GM12878_1))
hits3 <- findOverlaps(GRanges(PAX5_GM12878[-hits1@from,]),GRanges(CHIA_GM12878_1))
hits4 <- findOverlaps(GRanges(PAX5_GM12878[hits1@from,][unique(hits2@from),]),GRanges(hg19_promoter))
hits5 <- findOverlaps(GRanges(PAX5_GM12878[hits1@from,]),GRanges(CHIA_GM12878_1[-unique(idx1$V1),][-unique(idx4$V1),]))

write.table(PAX5_GM12878[hits1@from,][unique(hits5@from),][,c(1:3,5)],
            '/media/shaoqizhu/easystore/CTCF/GREAT/PAX5_unique_loop.bed',
            quote = F,col.names = F,row.names = F,sep = '\t')
hits6 <- findOverlaps(GRanges(PAX5_GM12878[hits1@from,]),GRanges(CHIA_GM12878_1))
hits7 <- findOverlaps(GRanges(PAX5_GM12878[hits1@from,]),GRanges(CHIA_GM12878_2))
length(unique(hits7@from[which(hits7@to %in% intersect(unique(hits6@to),unique(hits7@to)))]))
overlap <- unique(intersect(unique(hits6@to),unique(hits7@to)))
length(unique(hits6@from))
hits12 <- findOverlaps(GRanges(CHIA_GM12878_1[overlap,]),GRanges(hg19_promoter))
length(unique(findOverlaps(GRanges(CHIA_GM12878_2[overlap,][unique(hits12@from),]),GRanges(hg19_promoter))@from))

length(unique(findOverlaps(GRanges(CHIA_GM12878_2[overlap,]),GRanges(hg19_promoter))@from))


unique1 <- setdiff(unique(hits6@to),unique(hits7@to))
unique2 <- setdiff(unique(hits7@to),unique(hits6@to))
hits8 <- findOverlaps(GRanges(CHIA_GM12878_2[unique1,]),GRanges(CTCF_GM12878[-hits1@to,]))
hits9 <- findOverlaps(GRanges(CHIA_GM12878_2[unique1,][-hits8@from,]),GRanges(PAX5_GM12878))
length(unique(hits13@from))
hits13 <- findOverlaps(GRanges(CHIA_GM12878_1[unique1,]),GRanges(hg19_promoter))
length(unique(findOverlaps(GRanges(CHIA_GM12878_2[unique1,][unique(hits13@from),]),GRanges(hg19_promoter))@from))

hits10 <- findOverlaps(GRanges(CHIA_GM12878_1[unique2,]),GRanges(CTCF_GM12878[-hits1@to,]))
hits11 <- findOverlaps(GRanges(CHIA_GM12878_1[unique2,][-hits10@from,]),GRanges(PAX5_GM12878))
length(unique(hits14@from))
hits14 <- findOverlaps(GRanges(CHIA_GM12878_2[unique2,]),GRanges(hg19_promoter))
length(unique(findOverlaps(GRanges(CHIA_GM12878_1[unique2,][unique(hits14@from),]),GRanges(hg19_promoter))@from))

Unique <- as.numeric(rownames(CHIA_GM12878_1[-unique(idx1$V1),][-unique(idx4$V1),]))
hits6 <- findOverlaps(GRanges(PAX5_GM12878[hits1@from,]),GRanges(CHIA_GM12878_1[-Unique,]))
length(unique(hits6@from))

hits1 <- findOverlaps(GRanges(CEBPB_Hela),GRanges(CTCF_Hela))
hits2 <- findOverlaps(GRanges(CEBPB_Hela[hits1@from,]),GRanges(CHIA_Hela_1))
hits3 <- findOverlaps(GRanges(CEBPB_Hela[-hits1@from,]),GRanges(CHIA_Hela_1))
hits4 <- findOverlaps(GRanges(CEBPB_Hela[hits1@from,][unique(hits2@from),]),GRanges(hg19_promoter))
hits5 <- findOverlaps(GRanges(CEBPB_Hela[hits1@from,]),GRanges(CHIA_Hela_1[-unique(idx1$V2),][-unique(idx5$V1),]))
length(unique(hits5@from))
write.table(CEBPB_Hela[hits1@from,][unique(hits5@from),][,c(1:3,5)],
            '/media/shaoqizhu/easystore/CTCF/GREAT/CEBPB_unique_loop.bed',
            quote = F,col.names = F,row.names = F,sep = '\t')


hits1 <- findOverlaps(GRanges(GATA1_K562),GRanges(CTCF_K562))
hits2 <- findOverlaps(GRanges(GATA1_K562[hits1@from,]),GRanges(CHIA_K562_1))
hits3 <- findOverlaps(GRanges(GATA1_K562[-hits1@from,]),GRanges(CHIA_K562_1))
hits4 <- findOverlaps(GRanges(GATA1_K562[hits1@from,][unique(hits2@from),]),GRanges(hg19_promoter))
hits5 <- findOverlaps(GRanges(GATA1_K562[hits1@from,]),GRanges(CHIA_K562_1[-unique(idx3$V1),][-unique(idx6$V1),]))
length(unique(hits5@to))
write.table(GATA1_K562[hits1@from,][unique(hits5@from),][,c(1:3,5)],
            '/media/shaoqizhu/easystore/CTCF/GREAT/GATA1_unique_loop.bed',
            quote = F,col.names = F,row.names = F,sep = '\t')

length(unique(hits2@from))
length(unique(hits3@from))

length(unique(hits2@to))
length(unique(hits3@to))

findOverlaps(GRanges(GATA1_K562[hits1@from,][unique(hits2@from),]),GRanges(hg19_promoter))

findOverlaps(GRanges(TCF7L2_Hela),GRanges(CTCF_Hela_unique))
findOverlaps(GRanges(GATA1_K562),GRanges(CTCF_K562_prom))
findOverlaps(GRanges(BCL11A_GM12878),GRanges(CTCF_GM12878))
findOverlaps(GRanges(CEBPB_Hela),GRanges(CTCF_Hela))

CTCF_motif <- function(x){
motif_idx <- matchMotifs(human_pwms_v1, GRanges(x),genome = BSgenome.Hsapiens.UCSC.hg19)
match_num <- apply(motif_idx@assays@data$motifMatches,2,sum)
match_num[which(motif_idx@colData$name=='CTCF')]
which(motif_idx@assays@data$motifMatches[,which(motif_idx@colData$name=='CTCF')])
}
motif_idx_G <- CTCF_motif(CTCF_GM12878)
motif_idx_H <- CTCF_motif(CTCF_Hela)
motif_idx_K <- CTCF_motif(CTCF_K562)

CTCF_K562[CTCF_K562_motif_idx,]
Venn(CTCF_GM12878[motif_idx_G,],CTCF_Hela[motif_idx_H,],CTCF_K562[motif_idx_K,])
Venn(CTCF_GM12878,CTCF_Hela,CTCF_K562)

findOverlaps(GRanges(CTCF_GM12878_unique),GRanges(CHIA_GM12878_1))
findOverlaps(GRanges(CTCF_Hela_unique),GRanges(CHIA_Hela_1))
findOverlaps(GRanges(CTCF_K562_unique),GRanges(CHIA_K562_1))

hits1 <- findOverlaps(GRanges(CTCF_GM12878),GRanges(CTCF_Hela))
hits2 <- findOverlaps(GRanges(CTCF_GM12878[unique(hits1@from),]),GRanges(CTCF_K562))
hits3 <- findOverlaps(GRanges(CTCF_GM12878[unique(hits1@from),][unique(hits2@from),]),GRanges(CHIA_GM12878_1))
length(unique(hits3@to))


ESC_CHIA <- read.table('/media/shaoqizhu/easystore/CTCF/ESC/CHIA_PET_CTCF.bedpe')
ESC_CHIA_1 <- ESC_CHIA[,1:3]; colnames(ESC_CHIA_1) <- c('chr','start','end')
ESC_CHIA_2 <- ESC_CHIA[,4:6]; colnames(ESC_CHIA_2) <- c('chr','start','end')
ESC_OCT4 <- read.table('/media/shaoqizhu/easystore/CTCF/ESC/ESC_OCT4_GSM803438.bed')
ESC_SOX2 <- read.table('/media/shaoqizhu/easystore/CTCF/ESC/ESC_SOX2_GSM456570.bed')
ESC_CTCF <- read.table('/media/shaoqizhu/easystore/CTCF/ESC/ESC_CTCF_GSM733672.bed')
colnames(ESC_OCT4)[1:3] <- c('chr','start','end')
colnames(ESC_SOX2)[1:3] <- c('chr','start','end')
colnames(ESC_CTCF)[1:3] <- c('chr','start','end')

findOverlaps(GRanges(ESC_OCT4),GRanges(ESC_SOX2))
findOverlaps(GRanges(ESC_OCT4),GRanges(ESC_CTCF))
hits <- findOverlaps(GRanges(ESC_OCT4),GRanges(ESC_CTCF))
hits0 <- findOverlaps(GRanges(ESC_OCT4[hits@from,]),GRanges(ESC_SOX2))
findOverlaps(GRanges(ESC_OCT4[hits@from,][hits0@from,]),GRanges(hg38_promoter))
findOverlaps(GRanges(ESC_CTCF),GRanges(hg38_promoter))

hits1 <- findOverlaps(GRanges(ESC_CTCF),GRanges(ESC_CHIA_1))
hits2 <- findOverlaps(GRanges(ESC_CTCF),GRanges(ESC_CHIA_2))
loop_with_ctcf <- intersect(hits1@to,hits2@to)
unique(hits2@from)
findOverlaps(GRanges(ESC_CTCF[unique(hits1@from),]),GRanges(hg38_promoter))
length(unique(hits1@to))

hits3 <- findOverlaps(GRanges(ESC_OCT4[hits@from,]),GRanges(ESC_CTCF[hits2@from,]))
length(unique(hits3@from))
findOverlaps(GRanges(ESC_OCT4[hits@from,][unique(hits3@from),]),GRanges(hg38_promoter))

hits4 <- findOverlaps(GRanges(ESC_OCT4[-hits@from,]),GRanges(ESC_CHIA_1))
length(unique(hits4@to))
findOverlaps(GRanges(ESC_OCT4[-hits@from,][unique(hits4@from),]),GRanges(hg38_promoter))

ESC_CTCF_CHIA_1 <- cbind(ESC_CTCF[hits1@from,],from=hits1@from,ESC_CHIA_1[hits1@to,],to=hits1@to)
ESC_CTCF_CHIA_2 <- cbind(ESC_CTCF[hits2@from,],from=hits2@from,ESC_CHIA_2[hits2@to,],to=hits2@to)
findOverlaps(GRanges(ESC_CTCF[unique(ESC_CTCF_CHIA_1[ESC_CTCF_CHIA_1$to%in%loop_with_ctcf,]$from),]),GRanges(hg38_promoter))

ESC_CTCF_CHIA_2[ESC_CTCF_CHIA_2$to==loop_with_ctcf[1],]
ESC_CHIA[loop_with_ctcf[1],]

