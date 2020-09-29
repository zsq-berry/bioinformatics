library(openxlsx)
#sed -i 's/^[ \t]*//g' alignment.txt
PATH='/media/shaoqizhu/easystore/CD8-HP/CTCF/alignment/batch_2/'
alignment <- read.csv(paste(PATH,'alignment.txt',sep = ''),header = F,sep=' ')
alignment_bam_rmdup <- read.csv(paste(PATH,'bam.txt',sep = ''),header = F,sep=' ')
peaks <- read.csv(paste(PATH,'peaks.txt',sep = ''),header = F,sep=' ')
name <- read.csv(paste(PATH,'name.txt',sep = ''),header = F,sep=' ')
abstract <- as.data.frame(cbind(sample=as.character(name$V1),
                                total.reads=as.character(alignment$V1[seq(1,nrow(alignment),6)]),
                                mapping.rate=as.character(alignment$V1[seq(6,nrow(alignment),6)]),
                                `concordant.reads(rate)`=paste(as.character(alignment$V1[seq(4,nrow(alignment),6)]),as.character(alignment$V2[seq(4,nrow(alignment),6)]),sep=''),
                                bam.rmdup.reads=as.character(alignment_bam_rmdup$V1[seq(7,nrow(alignment_bam_rmdup),14)])
                                #peak.calling=as.character(peaks$V4)
                                ))


write.xlsx(abstract,paste(PATH,'abstract.xlsx',sep = ''))

alignment$V1[seq(4,nrow(alignment),24)]
abstract <- as.data.frame(cbind(sample=as.character(name$V1),
                                total.reads=as.character(alignment$V1[seq(4,nrow(alignment),24)]),
                                mapping.rate=as.character(alignment$V1[seq(22,nrow(alignment),24)]),
                                `concordant.reads(rate)`=paste(as.character(alignment$V1[seq(7,nrow(alignment),24)]),
                                                               as.character(alignment$V2[seq(7,nrow(alignment),24)]),sep=''),
                                bam.rmdup.reads=as.character(alignment_bam_rmdup$V1[seq(7,nrow(alignment_bam_rmdup),14)])
                                #peak.calling=as.character(peaks$V4)
))


rmdup <- as.data.frame(cbind(name=as.character(name$V1),
        bam.rmdup=as.character(alignment_bam_rmdup$V1[seq(7,nrow(alignment_bam_rmdup),14)])))


alignment <- read.csv('/media/shaoqizhu/easystore/LSD1/alignment/alignment-1.txt',header = F,sep=' ')
alignment_bam <- read.csv('/media/shaoqizhu/easystore/LSD1/alignment/alignment-bam.txt',header = F,sep=' ')
name <- read.csv('/media/shaoqizhu/easystore/LSD1/alignment/name.txt',header = F,sep=' ')
abstract <- as.data.frame(cbind(sample=as.character(name$V1),
                                total.reads=as.character(alignment$V1[seq(1,nrow(alignment),6)]),
                                mapping.rate=as.character(alignment$V1[seq(6,nrow(alignment),6)]),
                                `concordant.reads(rate)`=paste(as.character(alignment$V1[seq(4,nrow(alignment),6)]),as.character(alignment$V2[seq(4,nrow(alignment),6)]),sep=''),
                                bam.reads=as.character(alignment_bam$V1[seq(7,nrow(alignment_bam),14)])))
abstract
write.xlsx(abstract,'/media/shaoqizhu/easystore/LSD1/alignment/abstract.xlsx')
