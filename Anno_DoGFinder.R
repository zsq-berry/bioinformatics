library(stringr)
library(readr)
annotation <- read.csv('~/Downloads/ucsc/mm9_loci_annotation.bed',sep='\t',header = F)
name <- read.csv('~/Downloads/ucsc/mm9_name.bed',header = F)
ensembl <- read.csv('~/Downloads/ucsc/ensembl.bed',sep='\t')
refseq <- read.csv('~/Downloads/ucsc/refseq.bed',sep ='\t')
ucsc <- read.csv('~/Downloads/ucsc/ucsc.bed',sep ='\t')

ensembl_gene <- read.csv('~/Downloads/ucsc/ensembl_gene.bed',sep='\t')
refseq_gene <- read.csv('~/Downloads/ucsc/refseq_gene.bed',sep ='\t')
ucsc_gene <- read.csv('~/Downloads/ucsc/ucsc_gene.bed',sep ='\t')

ensembl_exonStarts <- str_split(ensembl_gene$exonStarts,",")
ensembl_exonEnds <- str_split(ensembl_gene$exonEnds,",")

ensembl$exonLen <- sapply(X = 1:length(ensembl_exonEnds), FUN = function(x){L = length(ensembl_exonEnds[[x]]);
  sum(as.numeric(ensembl_exonEnds[[x]][1:(L-1)]) - as.numeric(ensembl_exonStarts[[x]][1:(L-1)]))})

ucsc_exonStarts <- str_split(ucsc_gene$exonStarts,",")
ucsc_exonEnds <- str_split(ucsc_gene$exonEnds,",")

ucsc$exonLen <- sapply(X = 1:length(ucsc_exonEnds), FUN = function(x){L = length(ucsc_exonEnds[[x]]);
sum(as.numeric(ucsc_exonEnds[[x]][1:(L-1)]) - as.numeric(ucsc_exonStarts[[x]][1:(L-1)]))})

refseq_exonStarts <- str_split(refseq_gene$exonStarts,",")
refseq_exonEnds <- str_split(refseq_gene$exonEnds,",")

refseq$exonLen <- sapply(X = 1:length(refseq_exonEnds), FUN = function(x){L = length(refseq_exonEnds[[x]]);
sum(as.numeric(refseq_exonEnds[[x]][1:(L-1)]) - as.numeric(refseq_exonStarts[[x]][1:(L-1)]))})


annotation$ID <- str_split_fixed(annotation$V4,"&",2)[,1]
annotation <- annotation[order(annotation$ID),]

ensembl <- ensembl[order(ensembl$X.name),]
ucsc <- ucsc[order(ucsc$X.kgID),]
refseq <- refseq[order(refseq$name),]

ensembl_point <- which(annotation$ID=='ENSMUST00000175633')
refseq_point <- which(annotation$ID=='NR_153101')

nrow(ensembl[ensembl$X.name %in% annotation$ID,])
setdiff(annotation$ID[1:ensembl_point],ensembl$X.name)
setdiff(annotation$ID[(ensembl_point+1):refseq_point],refseq$name)
setdiff(annotation$ID[(refseq_point+1):37504],ucsc$X.kgID)

annotation <- annotation[annotation$ID!='NR_105819',]
dupname <- names(table(annotation$ID)[table(annotation$ID)>1])
for (i in dupname){
  annotation <- annotation[annotation$ID!=i,]
}

anno_ref <- as.data.frame(refseq[refseq$name%in%annotation$ID[(ensembl_point+1):refseq_point],c("name","name2")])
dupname1 <- names(table(as.character(anno_ref$name))[table(as.character(anno_ref$name))>1])
for (i in dupname1){
  annotation <- annotation[annotation$ID!=i,]
}

annotation$name <- NA
annotation$exonlen <- NA

ensembl_point <- which(annotation$ID=='ENSMUST00000175633')
refseq_point <- which(annotation$ID=='NR_153101')

annotation$name[1:ensembl_point] <- as.character(ensembl$value[ensembl$X.name%in%annotation$ID])
annotation$name[(ensembl_point+1):refseq_point] <- unique(as.character(refseq$name2[refseq$name%in%annotation$ID[(ensembl_point+1):refseq_point]]))
annotation$name[(refseq_point+1):nrow(annotation)] <- as.character(ucsc$geneSymbol[ucsc$X.kgID%in%annotation$ID])

annotation$exonlen[1:ensembl_point] <- as.numeric(ensembl$exonLen[ensembl$X.name%in%annotation$ID])
annotation$exonlen[(ensembl_point+1):refseq_point] <- as.numeric(refseq$exonLen[refseq$name%in%annotation$ID[(ensembl_point+1):refseq_point]])
annotation$exonlen[(refseq_point+1):nrow(annotation)] <- as.numeric(ucsc$exonLen[ucsc$X.kgID%in%annotation$ID])


anno <- annotation[,c("V1","V2","V3","V6","name","exonlen")]

if(length(table(table(anno$name)))>1){
  anno <- anno[order(anno$name),]
  anno.rep <- as.data.frame(table(anno$name)[table(anno$name)>1])
  idx=as.numeric(); for(i in 1:nrow(anno.rep)){idx = c(idx,0:(anno.rep$Freq[i]-1))}
  rep=anno$name %in% anno.rep$Var1
  anno$name[rep][idx!=0]=paste(anno$name[rep][idx!=0],idx[idx!=0],sep='-')
}
anno$name <- gsub("_","-",anno$name)
write_delim(anno,paste('/media/shaoqizhu/easystore/read-through/annotation/','/','geneAnno.bed',sep=''),delim = "\t",col_names=F)
