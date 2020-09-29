ctcf_ctrl<-list()
ctcf_auxin1d<-list()
ctcf_auxin2d<-list()
ctcf_auxin4d<-list()
ctcf_washoff2d<-list()
for(i in 1:3){
  ctcf_ctrl[[i]] <- read.csv(paste('/media/shaoqizhu/easystore/read-through/ctcf/ctcf_ctrl_',i,'_readthrough.csv',sep=''))
  ctcf_auxin1d[[i]] <- read.csv(paste('/media/shaoqizhu/easystore/read-through/ctcf/ctcf_auxin1d_',i,'_readthrough.csv',sep=''))
  ctcf_auxin2d[[i]] <- read.csv(paste('/media/shaoqizhu/easystore/read-through/ctcf/ctcf_auxin2d_',i,'_readthrough.csv',sep=''))
  ctcf_auxin4d[[i]] <- read.csv(paste('/media/shaoqizhu/easystore/read-through/ctcf/ctcf_auxin4d_',i,'_readthrough.csv',sep=''))
  ctcf_washoff2d[[i]] <- read.csv(paste('/media/shaoqizhu/easystore/read-through/ctcf/ctcf_washoff2d_',i,'_readthrough.csv',sep=''))
}

i=2; sample_1 <- ctcf_ctrl[[i]]; sample_2 <- ctcf_washoff2d[[i]]
rt_overlap <- intersect(sample_1$gene,sample_2$gene)
hist(log((sample_2$rpkm/sample_2$gene.rpkm)[sample_2$gene%in%rt_overlap])-
     log((sample_1$rpkm/sample_1$gene.rpkm)[sample_1$gene%in%rt_overlap]),
     main=paste('readthrough rpkm','rep',i),xlab='log(auxin1d/ctrl)')


i=2; sample_1 <- ctcf_ctrl[[i]]; sample_2 <- ctcf_auxin1d[[i]]
rt_overlap <- intersect(sample_1$gene,sample_2$gene)
plot(log((sample_1$len)[sample_1$gene%in%rt_overlap]),
       log((sample_2$len)[sample_2$gene%in%rt_overlap]),main='replicate 2 rt length',xlab='ctrl',ylab='ctcf_auxin2d')
hist(log((sample_2$len)[sample_2$gene%in%rt_overlap])-
       log((sample_1$len)[sample_1$gene%in%rt_overlap]),
     main=paste('readthrough rpkm','rep',i),xlab='log(auxin4d/ctrl)')
which(log((sample_2$len)[sample_2$gene%in%rt_overlap])-
        log((sample_1$len)[sample_1$gene%in%rt_overlap])>0)

sample_2$gene[sample_2$gene%in%rt_overlap][which(sample_2$len[sample_2$gene%in%rt_overlap]-
                                                   sample_1$len[sample_1$gene%in%rt_overlap]>(5))]
sample_1$gene[sample_1$gene%in%rt_overlap][which(sample_2$len[sample_2$gene%in%rt_overlap]-
                                                   sample_1$len[sample_1$gene%in%rt_overlap]>(5))]

which(rt_overlap=="Mrps5")

geneAnno <- read.table('/media/shaoqizhu/easystore/read-through/annotation/geneAnno.bed')
geneAnno$V6[geneAnno$V5%in%which(log((sample_2$len)[sample_2$gene%in%rt_overlap])-
                                 log((sample_1$len)[sample_1$gene%in%rt_overlap])<0)]
sample_1[-which(sample_1$gene%in%rt_overlap),][order(sample_1[-which(sample_1$gene%in%rt_overlap),]$rpkm,decreasing = T),]

mm10 <- read.table('/media/shaoqizhu/easystore/Annotation/mm10_2015.bed')
mm10$v6 <- 1:nrow(mm10)
mm10$V7 <- mm10$V3-mm10$V2
mm10 <- mm10[,c(1,2,3,5,6,4,7)]
write.table(mm10,'/media/shaoqizhu/easystore/read-through/annotation/geneAnno.bed',sep = '\t',col.names = F,row.names = F,quote = F)


ctrl <- read.csv(paste('/media/shaoqizhu/easystore/read-through/ctcf/Ctrl-1_readthrough.csv',sep=''))
HS_1 <- read.csv(paste('/media/shaoqizhu/easystore/read-through/ctcf/HS-1_readthrough.csv',sep=''))
H2O2_1 <- read.csv(paste('/media/shaoqizhu/easystore/read-through/ctcf/H2O2-1_readthrough.csv',sep=''))
KCl_1 <- read.csv(paste('/media/shaoqizhu/easystore/read-through/ctcf/KCl-1_readthrough.csv',sep=''))

sample_1 <- ctrl; sample_2 <- HS_1
rt_overlap <- intersect(sample_1$gene,sample_2$gene)
hist(log((sample_2$len)[sample_2$gene%in%rt_overlap])-
       log((sample_1$len)[sample_1$gene%in%rt_overlap]),
     main=paste('readthrough rpkm','rep',i),xlab='log(auxin1d/ctrl)')
plot(log((sample_2$len)[sample_2$gene%in%rt_overlap]),
     log((sample_1$len)[sample_1$gene%in%rt_overlap]))

ctcf2_auxin0h<-list()
ctcf2_auxin1d<-list()
ctcf2_auxin2d<-list()
ctcf2_auxin8d<-list()
for(i in 1:2){
ctcf2_auxin0h[[i]] <- read.csv(paste('/media/shaoqizhu/easystore/read-through/ctcf2/auxin_0h_rep',i,'_readthrough.csv',sep=''))
ctcf2_auxin1d[[i]] <- read.csv(paste('/media/shaoqizhu/easystore/read-through/ctcf2/auxin_24h_rep',i,'_readthrough.csv',sep=''))
ctcf2_auxin2d[[i]] <- read.csv(paste('/media/shaoqizhu/easystore/read-through/ctcf2/auxin_48h_rep',i,'_readthrough.csv',sep=''))
ctcf2_auxin8d[[i]] <- read.csv(paste('/media/shaoqizhu/easystore/read-through/ctcf2/auxin_96h_rep',i,'_readthrough.csv',sep=''))
}
  
i=1; sample_1 <- ctcf2_auxin0h[[i]]; sample_2 <- ctcf_auxin2d[[i]]
rt_overlap <- intersect(sample_1$gene,sample_2$gene)
plot(log((sample_1$len)[sample_1$gene%in%rt_overlap]),
     log((sample_2$len)[sample_2$gene%in%rt_overlap]),main='replicate 2 rt length',xlab='ctrl',ylab='washoff2d')
hist(log((sample_2$len)[sample_2$gene%in%rt_overlap])-
       log((sample_1$len)[sample_1$gene%in%rt_overlap]),
     main=paste('readthrough rpkm','rep',i),xlab='log(auxin4d/ctrl)')
