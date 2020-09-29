library(readr)
library(rtracklayer)
library(msm)
library(expm)
library(stringr)
H2O2_mtx <- read.csv('/media/shaoqizhu/easystore/read-through/bigwig/mm9/H2O2-1_plus_all.mat',header = F,sep = '\t')
H2O2_mtx$mean <- apply(H2O2_mtx[,7:26],1,mean)


GM23338_mtx <- read.csv('/media/shaoqizhu/easystore/read-through/bigwig/hg19/fibroblast_dermis_minus_ENCFF786STS.mat',header = F,sep = '\t')
GM23338_mtx$mean <- apply(GM23338_mtx[,7:16],1,mean)


CD4_mtx <- read.csv('/media/shaoqizhu/easystore/read-through/bigwig/hg19/CD4/CD4_plus.mat',header = F,sep = '\t')
CD4_mtx$mean <- apply(CD4_mtx[,7:26],1,mean)

intergenic_plus <- read.table('/media/shaoqizhu/easystore/read-through/annotation/intergenic_plus.bed')
colnames(intergenic_plus)[1:3] <- c('chr','start','end')
Ctrl_plus <- import.bw('/media/shaoqizhu/easystore/read-through/bigwig/mm9/Ctrl-1_plus.bw')
H2O2_plus <- import.bw('/media/shaoqizhu/easystore/read-through/bigwig/mm9/H2O2-1_plus.bw')

gene <- which(intergenic_plus$start==61821916)
hits <- findOverlaps(GRanges(intergenic_plus),Ctrl_plus,minoverlap = 100)

as.data.frame(Ctrl_plus[hits@to])

plot(Ctrl_plus$score[hits@to[hits@from==6]])

intergenic_coverage <- as.data.frame(cbind(hits@from,Ctrl_plus$score[hits@to],Ctrl_plus@ranges@width[hits@to]/100))
sum(intergenic_coverage$V3)


aggregate(intergenic_coverage$V3,list(intergenic_coverage$V1),sum)[2]

rt_stats = function(rc){
  rt.stats <- aggregate(rc$V2,list(rc$V1),sum)
  rt.stats[3] <- aggregate(rc$V3,list(rc$V1),sum)[2]
  rt.stats[2] <- rt.stats[2]/rt.stats[3]
  rt.stats[4] <- aggregate(rc$V2,list(rc$V1),FUN = function(x){sum(x>=1)/length(x)})[2]
  rt.stats[5] <- aggregate(rc$V2,list(rc$V1),FUN = function(x){if(length(x)>1){max(abs(x[1:length(x)-1]-x[2:length(x)]))}else{0}})[2]
  rt.stats[6] <- aggregate(rc$V2,list(rc$V1),FUN = function(x){if(length(x)>=5){sum(x[1:5]>=1)}else{sum(x>=1)}})[2]
  names(rt.stats) = c('name','mean','len','count_thred','max_change','tes_cover')
  rt.stats
}
intergenic_stat <- rt_stats(intergenic_coverage)

gene_train = function(rt.stats){
  genes <- rt.stats[rt.stats$mean>0.5 & rt.stats$mean<50 
                    & rt.stats$len>=2 & rt.stats$max_change>0 & rt.stats$max_change<50,]
  genes <- genes[order(genes$mean,decreasing = T),]
  if (nrow(genes)>100) {intergenic_train <- intergenic_coverage[intergenic_coverage$V1 %in% genes$name[sample(1:nrow(genes),100)],]}
  else {intergenic_train <- intergenic_coverage[intergenic_coverage$V1 %in% genes$name,]}
}

intergenic_train <- gene_train(intergenic_stat)
sum(intergenic_train$V3)

intergenic_extended <- as.numeric();intergenic_idx <- as.numeric()
for (i in 1:nrow(intergenic_train)){
  intergenic_extended <- c(intergenic_extended,rep(intergenic_train$V2[i],intergenic_train$V3[i]))
  intergenic_idx <- c(intergenic_idx,rep(intergenic_train$V1[i],intergenic_train$V3[i]))
}

intergenic_train <- as.data.frame(cbind(V1=intergenic_idx,V2=intergenic_extended))

unique(intergenic_train$V1)

states = intergenic_train$V2
timeq=1:length(states);m = rbind(c(0.9,0.1),c(0,1))
hmm.model <- list(hmmNorm(mean = mean(intergenic_train$V2),sd =sd(intergenic_train$V2)),hmmNorm(mean=0.1, sd=1));
fitted.msm <- msm(states~timeq,subject = intergenic_train$V1,
                  qmatrix = logm(m), hmodel = hmm.model)

print(fitted.msm$paramdata$estimates.t)
print(pmatrix.msm(fitted.msm))
print(logLik.msm(fitted.msm))

integenic_vit <- viterbi.msm(fitted.msm)[,4]

rt <- unique(intergenic_train$V1)[28]
intergenic <- which(intergenic_train$V1==rt)
plot(intergenic_train$V2[intergenic])
plot(2-integenic_vit[intergenic])
sum(2-integenic_vit[intergenic])

intergenic_plus[rt,]

start <- Ctrl_plus@ranges@start[hits@to[hits@from==rt]][1]
chr <- Ctrl_plus@seqnames[hits@to[hits@from==rt]][1]
as.character(chr)
cbind(start,start+sum(2-integenic_vit[intergenic])*100)




# library(MASS)
# para <- matrix(data = 0,nrow = 100,ncol = 2)
# for (i in 1:100){
#   intergenic <- intergenic_train$V2[intergenic_train$V1==unique(intergenic_train$V1)[i]]
#   fit_region <- ifelse(length(intergenic)>=6,floor(length(intergenic)/6),length(intergenic))
#   para[i,] <- fitdistr(intergenic[1:fit_region],'normal')$estimate
# }

library(gridExtra)
library(cowplot)
library(ggplot2)

gene_name <- unique(intergenic_vit$V1)
plot_gene(intergenic_vit,5)

plot_gene <- function(rt,N){
  g<-list()
  for (i in 1:10){
    Data = rt[rt$V1==gene_name[i+10*N],]; Data$V3 = 1:nrow(Data)
    g[[i]] <- ggplot() + geom_line(data=Data,aes(x=V3,y=V2),colour="darkgrey") +
      geom_vline(xintercept=sum(Data$vit),linetype=4,colour="blue") +
      labs(x = "bin's position",y="reads count",title= gene_name[i+10*N]) +
      theme(plot.title=element_text(size=12,vjust=0.5,hjust=0.5))}
  plot_grid(plotlist = g, ncol = 5)}


library(GenomicFeatures)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
mm9_ucsc <- TxDb.Mmusculus.UCSC.mm9.knownGene
transcriptsBy(mm9_ucsc,by='gene')
columns(mm9_ucsc)
columns(org.Mm.eg.db)
ucsc_id <- keys(mm9_ucsc,'TXNAME')
df = select(mm9_ucsc, ucsc_id, "GENEID",'TXNAME')
ucsc <- select(org.Mm.eg.db, df$GENEID, c("SYMBOL"))
ucsc$txname <- df$TXNAME

columns(org.Mm.eg.db)

ensembl <- select(org.Mm.eg.db, df$GENEID, c("ENSEMBLTRANS","SYMBOL"))
refseq <- select(org.Mm.eg.db, df$GENEID, c("REFSEQ","SYMBOL"))

name <- c("Ctrl-1","Ctrl-2","H2O2-1","H2O2-2","HS-1","HS-2","KCl-1","KCl-2")
dog_finder <- list()
for (i in 1:8){
  dog_finder[[i]] <- read.table(paste('/media/shaoqizhu/easystore/read-through/sam/dogFinder/dogFinder/',name[i],'/Final_Dog_annotation.bed',sep = ''),stringsAsFactors = F)
}
names(dog_finder) <- gsub('-','_',name)
txname_list <- str_split(dog_finder$Ctrl_1$V4,"&")
txname <- as.data.frame(cbind(unlist(txname_list),rep(1:length(txname_list),times=lengths(txname_list))))
setdiff(txname$V1[grep('EN',txname$V1)],ensembl$ENSEMBLTRANS)

require("biomaRt")
listMarts()
listDatasets(mart)
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("mmusculus_gene_ensembl", mart)

grep('ensembl_transcript_id',listFilters(mart)[1])
listFilters(mart)

listAttr<- listAttributes(mart)
listAttr$name[grep('uc',listAttr$name)]

Ensembl <- getBM(
  mart=mart,
  attributes=c("ensembl_transcript_id","external_gene_name"),
  filter="ensembl_transcript_id",
  values=txname$V1[grep('EN',txname$V1)])
txname$symble[grep('EN',txname$V1)] <- Ensembl$external_gene_name

refseq_NM <- getBM(
  mart=mart,
  attributes=c("refseq_mrna","external_gene_name"),
  filter="refseq_mrna",
  values=txname$V1[grep('NM',txname$V1)])
txname$symble[grep('NM',txname$V1)] <- refseq_NM$external_gene_name

refseq_NR <- getBM(
  mart=mart,
  attributes=c("refseq_ncrna","external_gene_name"),
  filter="refseq_ncrna",
  values=txname$V1[grep('NR',txname$V1)])
txname$symble[grep('NR',txname$V1)] <- refseq_NR$external_gene_name

txname <- txname[order(txname$V1),]
ucsc <- ucsc[order(ucsc$txname),]
txname$symble[grep('uc',txname$V1)] <- ucsc[ucsc$txname %in% txname$V1[grep('uc',txname$V1)],2]


name <- c("Ctrl-1","Ctrl-2","H2O2-1","H2O2-2","HS-1","HS-2","KCl-1","KCl-2")
RTtools <- list()
for (i in 1:8){
  RTtools[[i]] <- read.csv(paste('/media/shaoqizhu/easystore/read-through/sam/dogFinder/RTtools/',name[i],'_readthrough.csv',sep = ''),stringsAsFactors = F)
}
names(RTtools) <- gsub('-','_',name)


library(readxl)
dogFinder <- read_excel('/media/shaoqizhu/easystore/read-through/dogFinder.xls')

setdiff(dogFinder$`gene name`[dogFinder$UN...8!='NaN'],
        RTtools$Ctrl_1$gene_name[RTtools$Ctrl_1$rt_length>=20])
gene <- setdiff(RTtools$Ctrl_1$gene_name[RTtools$Ctrl_1$rt_length>=20],
        dogFinder$`gene name`[dogFinder$UN...8!='NaN'])

setdiff(dogFinder$`gene name`[dogFinder$H2O2...9!='NaN'],
        RTtools$H2O2_1$gene_name[RTtools$H2O2_1$rt_length>=20])

dogFinder[dogFinder$UN...8!='NaN']

RTtools$Ctrl_1[RTtools$Ctrl_1$gene_name%in% gene,][200:300,]



mm9_genome <- read.table('/media/shaoqizhu/easystore/read-through/Rcode/gtf/mm9.genome')

write.table(cbind(as.character(mm9_genome$V1),0,mm9_genome$V2,as.character(mm9_genome$V1)),
                  '/media/shaoqizhu/easystore/read-through/Rcode/gtf/mm9_genome.bed',
                  quote = F,row.names = F,col.names = F,sep = '\t')
#cat mm9_genome.bed | sort -k1,1 -k2,2n | bedtools makewindows -b stdin -w 200 -i srcwinnum | gzip > mm9_windows.bed.gz
#cat Ctrl-1_sorted_pe.bed | grep -w + | bedtools intersect -a ../../Rcode/gtf/mm9_windows.bed.gz -b stdin -wao -sorted | awk '$11>0 {$11=$11/($7-$6)} {print $1"\t"$2"\t"$3"\t"$11}'| gzip > Ctrl-1_coverage_plus.bed.gz








intergenic <- read.table('/media/shaoqizhu/easystore/read-through/Rcode/intergenic.bed')
intergenic <- intergenic[order(intergenic[,5]),]

gene_coverage <- read.table('/media/shaoqizhu/easystore/read-through/sam/dogFinder/RTtools/Ctrl-1_gene_coverage.bed')
gene_coverage <- gene_coverage[order(gene_coverage[,4]),]

cat("  Loading reads coverage of intergenic regions",sep="\n")
rt.coverage.plus = read.table('/media/shaoqizhu/easystore/read-through/sam/dogFinder/RTtools/Ctrl-1_coverage_plus.bed.gz')
rt.coverage.minus = read.table('/media/shaoqizhu/easystore/read-through/sam/dogFinder/RTtools/Ctrl-1_coverage_minus.bed.gz')

rt.coverage = rbind(rt.coverage.plus,rt.coverage.minus[nrow(rt.coverage.minus):1,])

#=====================Get coverage statistics for each gene's readthrough
#stats is short for statistics

rt_stats = function(rc){
  rt.stats <- aggregate(rc[,2],list(rc[,1]),mean)
  rt.stats[3] <- aggregate(rc[,2],list(rc[,1]),length)[2]
  rt.stats[4] <- aggregate(rc[,2],list(rc[,1]),FUN = function(x){sum(x>=1)/length(x)})[2]
  rt.stats[5] <- aggregate(rc[,2],list(rc[,1]),FUN = function(x){if(length(x)>1){max(abs(x[1:length(x)-1]-x[2:length(x)]))}else{0}})[2]
  rt.stats[6] <- aggregate(rc[,2],list(rc[,1]),FUN = function(x){if(length(x)>=5){sum(x[1:5]>=1)}else{sum(x>=1)}})[2]
  names(rt.stats) = c('name','mean','len','count_thred','max_change','tes_cover')
  rt.stats
}

#=====================Select proper readthrough genes for HMM training and all valid readthrough genes for viterbi fitting

gene_vit = function(rt.stats){
  # gene.valid <- rt.stats[rt.stats$mean>0.5 & rt.stats$mean<50 
  #                        & rt.stats$len>=2 & rt.stats$max_change>0
  #                        & rt.stats$max_change<30,"name"]
  gene.valid <- rt.stats[rt.stats$len>=5 & rt.stats$tes_cover>=2,"name"]
  rt.vit <- rt.coverage[rt.coverage[,1] %in% gene.valid,]
}

intergenic_stat <- rt_stats(rt.coverage)

intergenic_vit <- gene_vit(intergenic_stat)

hmm_fit <- function(states){
  tryCatch({
  timeq=1:length(states);m = rbind(c(0.9,0.1),c(0,1))
  hmm.model <- list(hmmNorm(mean = 5,sd = 1),hmmNorm(mean=0.4, sd=0.4));
  fitted.msm <- suppressWarnings(msm(states~timeq,qmatrix = logm(m), hmodel = hmm.model, fixedpars = c(4,5)))
  sum(2-viterbi.msm(fitted.msm)[,4])
  }, error=function(e){})
}

idx <- which(intergenic_vit$V1 %in% unique(intergenic_vit$V1)[1:2000])

{s=Sys.time()
rt_len <- aggregate(intergenic_vit$V2[idx],list(intergenic_vit$V1[idx]),hmm_fit)
e=Sys.time()
print(e-s)}




# =====================Training using msm package
suppressMessages(library(progress))
valid <- unique(intergenic_vit[,1])
pb <- progress_bar$new(format = "  HMM fitting [:bar] :percent eta: :eta",
                       total = length(valid), width= 70, clear = F)

intergenic_vit$vit <- 0
for (i in 1:length(valid)){
  tryCatch({
    gene <- which(intergenic_vit[,1]==valid[i])
    states = intergenic_vit[,2][gene]
    timeq=1:length(states);m = rbind(c(0.9,0.1),c(0,1))
    hmm.model <- list(hmmNorm(mean = mean(states),sd = sd(states)),hmmNorm(mean=0.4, sd=0.4));
    fitted.msm <- suppressWarnings(msm(states~timeq,qmatrix = logm(m), hmodel = hmm.model, fixedpars = c(4,5)))
    intergenic_vit$vit[gene] <- 2-viterbi.msm(fitted.msm)[,4]
  }, error=function(e){})
  pb$tick()
  Sys.sleep(1 / 100)
}
# pbsapply()
aggregate()




