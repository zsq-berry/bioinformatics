library(rtracklayer)
library(readr)

mm10 <- import('~/Documents/annotation/ucsc.mm10.annotation.gtf')

mm10 <- as.data.frame(mm10)
mm10_gene <- cbind(seqnames <- as.character(unlist(aggregate(as.character(mm10$seqnames),list(mm10$gene_name),FUN = function(x){x[1]})[,2])),
                   start <- aggregate(mm10$start,list(mm10$gene_name),min)[,2],
                   end <- aggregate(mm10$end,list(mm10$gene_name),max)[,2],
                   strand <- as.character(unlist(aggregate(as.character(mm10$strand),list(mm10$gene_name),FUN = function(x){x[1]})[,2])),
                   gene <- aggregate(mm10$start,list(mm10$gene_name),min)[,1])



refseq <- import('~/Documents/annotation/GCF_000001635.26_GRCm38.p6_genomic.gff.gz')
refseq <- as.data.frame(refseq)

gencode <- import('~/Documents/annotation/gencode.vM23.annotation.gtf.gz')
gencode <- as.data.frame(gencode)

promoter <- geneAnno
promoter$start <- geneAnno$start-1000
promoter$end <- geneAnno$start+1000

write.table(promoter,paste('/media/shaoqizhu/easystore/hdac/promoter.bed',sep=''),
            quote = F,sep = '\t',row.names = F,col.names = F)

peak_gene <- read.csv('~/Desktop/peaks/dKO_D0-WT_D0_positive.txt',sep = '\t',header = F)
scRNA_gene <- read.csv("~/Desktop/peaks/dKO_D0-WT_D0.csv")
length(scRNA_gene$X[scRNA_gene$avg_logFC>0])
gene <- intersect(peak_gene$V1,scRNA_gene$X[scRNA_gene$avg_logFC>0])
write.csv(gene,'~/Desktop/gene.csv')


great_negative <- read.csv('~/Desktop/peaks/dKO_D1-WT_D1_negative.txt',sep = '\t',header = F)
great_positive <- read.csv('~/Desktop/peaks/dKO_D1-WT_D1_positive.txt',sep = '\t',header = F)
updown50_gene <- read.csv('~/Desktop/summit/50kb/dKO_D1-WT_D1_gene.bed',sep = '\t',header = F)
Aa = length(intersect(great_negative$V1,updown50_gene$V6[updown50_gene$V5<0]))
Ab = length(intersect(great_negative$V1,updown50_gene$V6[updown50_gene$V5>0]))
Ba = length(intersect(great_positive$V1,updown50_gene$V6[updown50_gene$V5<0]))
Bb = length(intersect(great_positive$V1,updown50_gene$V6[updown50_gene$V5>0]))

A = length(great_negative$V1)
B = length(great_positive$V1)
a = length(updown50_gene$V6[updown50_gene$V5<0])
b = length(updown50_gene$V6[updown50_gene$V5>0])
as.data.frame(cbind(c(Aa,Ab,A),c(Ba,Bb,B),c(a,b,NA)))


updown50_gene <- read.csv('~/Desktop/summit/50kb/WT_D0-WT_D1_gene.bed',sep = '\t',header = F)
scRNA_gene <- read.csv("~/Desktop/peaks/WT_D0-WT_D1.csv")
Aa = (intersect(scRNA_gene$X[scRNA_gene$avg_logFC<0],updown50_gene$V6[updown50_gene$V5<0]))
Ab = (intersect(scRNA_gene$X[scRNA_gene$avg_logFC<0],updown50_gene$V6[updown50_gene$V5>0]))
Ba = (intersect(scRNA_gene$X[scRNA_gene$avg_logFC>0],updown50_gene$V6[updown50_gene$V5<0]))
Bb = (intersect(scRNA_gene$X[scRNA_gene$avg_logFC>0],updown50_gene$V6[updown50_gene$V5>0]))

A = length(scRNA_gene$X[scRNA_gene$avg_logFC<0])
B = length(scRNA_gene$X[scRNA_gene$avg_logFC>0])
a = length(updown50_gene$V6[updown50_gene$V5<0])
b = length(updown50_gene$V6[updown50_gene$V5>0])
as.data.frame(cbind(c(Aa,Ab,A),c(Ba,Bb,B),c(a,b,NA)))
intersect(Bb,TF$Symbol)

great_negative <- read.csv('/media/shaoqizhu/easystore/hdac/bampe/diff_peaks/great/dKO_D1_mt_WT_D1_great.bed',sep = '\t',header = F)
great_positive <- read.csv('/media/shaoqizhu/easystore/hdac/bampe/diff_peaks/great/dKO_D1_lt_WT_D1_great.bed',sep = '\t',header = F)
scRNA_gene <- read.csv("/media/shaoqizhu/easystore/hdac/bampe/diff_peaks/great/dKO_D1-WT_D1.csv")
Aa = length(intersect(great_negative$V1,scRNA_gene$X[scRNA_gene$avg_logFC<0]))
Ab = length(intersect(great_negative$V1,scRNA_gene$X[scRNA_gene$avg_logFC>0]))
Ba = length(intersect(great_positive$V1,scRNA_gene$X[scRNA_gene$avg_logFC<0]))
Bb = length(intersect(great_positive$V1,scRNA_gene$X[scRNA_gene$avg_logFC>0]))

A = length(great_negative$V1)
B = length(great_positive$V1)
a = length(scRNA_gene$X[scRNA_gene$avg_logFC<0])
b = length(scRNA_gene$X[scRNA_gene$avg_logFC>0])

as.data.frame(cbind(c(Aa,Ab,A),c(Ba,Bb,B),c(a,b,NA)))

TF <- read.csv('/Volumes/ShaoqiZhu/GenomicData/AnimalTFDB/Mus_musculus_TF.txt',sep = '\t')
intersect(Ab,TF$gene)
enrichr <- read.csv('~/Desktop/peaks/enrichr_WT_D1>dKO_D1.txt',sep='\t')
great_myc<- read.csv('~/Desktop/peaks/Myc_CD8_great.txt',sep = '\t',header = F)
intersect(toupper(great_myc$V1),str_split(enrichr$Genes[2],';')[[1]])
intersect(great_myc$V1,scRNA_gene$X[scRNA_gene$avg_logFC<0])
intersect(toupper(scRNA_gene$X[scRNA_gene$avg_logFC<0]),str_split(enrichr$Genes[2],';')[[1]])

k27ac_D0 <- read.csv('/Volumes/ShaoqiZhu/hdac1/bam/f_D0_K27ac_sorted_rmdup-W200-G600.scoreisland',header = F,sep = '\t')
k27ac_na <- read.csv('/Volumes/ShaoqiZhu/hdac1/bam/f_na_K27ac_sorted_rmdup-W200-G600.scoreisland',header = F,sep = '\t')
plot(sort(k27ac_D0$V4),main='K27ac_D0 island score ranking',xlab='index',ylab='score')


great_dKO_D0 <- read.csv('~/Desktop/peaks/dKO_D0-na_WT_D0_great.txt',sep = '\t')
great_WT_D0 <- read.csv('~/Desktop/peaks/WT_D0-na_dKO_D0_great.txt',sep = '\t')
scRNA_gene <- read.csv("~/Desktop/peaks/dKO_D0-WT_D0.csv")
Aa = length(intersect(great_dKO_D0[,1],scRNA_gene$X[scRNA_gene$avg_logFC<0]))
Ab = length(intersect(great_dKO_D0[,1],scRNA_gene$X[scRNA_gene$avg_logFC>0]))
Ba = length(intersect(great_WT_D0[,1],scRNA_gene$X[scRNA_gene$avg_logFC<0]))
Bb = length(intersect(great_WT_D0[,1],scRNA_gene$X[scRNA_gene$avg_logFC>0]))
A = length(great_dKO_D0[,1])
B = length(great_WT_D0[,1])
a = length(scRNA_gene$X[scRNA_gene$avg_logFC<0])
b = length(scRNA_gene$X[scRNA_gene$avg_logFC>0])
as.data.frame(cbind(c(Aa,Ab,A),c(Ba,Bb,B),c(a,b,NA)))

dKO_D0_gene <- intersect(rownames(HdacD01@assays$RNA@data),great_dKO_D0[,1])
WT_D0_gene <- intersect(rownames(HdacD01@assays$RNA@data),great_WT_D0[,1])
dKO_expr <- apply(HdacD01@assays$RNA@data[WT_D0_gene,HdacD0_Dko.data@Dimnames[[2]]],1,mean)
WT_expr <- apply(HdacD01@assays$RNA@data[WT_D0_gene,HdacD0.data@Dimnames[[2]]],1,mean)
plot(dKO_expr,WT_expr,main=paste('dKO_D0_gene',round(cor(dKO_expr,WT_expr),2)))
round(cor(dKO_expr,WT_expr),2)

expr_diff <- log2((dKO_expr+1)/(WT_expr+1))

which(expr_diff>(0.15))



