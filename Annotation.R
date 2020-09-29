library(rtracklayer)

gtf <- import('/media/shaoqizhu/easystore/Annotation/mm9_2015.gtf')
gtf <- as.data.frame(gtf)

gene.name = sort(unique(gtf$gene_name))
gtf_reduced <- as.data.frame(matrix(NA,ncol = 6,nrow = 0))
for (x in 1){
    gene = gtf[gtf$gene_name==gene.name[x],]
    transcript_chr = aggregate(gene$seqnames,list(gene$transcript_id),unique)[,2]
    transcript_start = aggregate(gene$start,list(gene$transcript_id),min)[,2]
    transcript_end = aggregate(gene$end,list(gene$transcript_id),max)[,2]
    transcript_strand = as.character(aggregate(gene$strand,list(gene$transcript_id),unique)[,2])
    transcript_ranges <- as.data.frame(reduce(GRanges(transcript_chr,IRanges(transcript_start,transcript_end),transcript_strand)))
    gene_exon = gene[gene$type=='exon',]; exonLen <- 0
    for (i in 1:nrow(transcript_ranges)){
      start_id <- which(gene_exon$start == transcript_ranges$start[i])
      end_id <- which(gene_exon$end == transcript_ranges$end[i])
      exonLen[i] <- sum(as.data.frame(reduce(GRanges(gene_exon[start_id[1]:end_id[length(end_id)],])))$width)
    }
    gtf_reduced = rbind(gtf_reduced,cbind(transcript_ranges,gene_name=rep(gene.name[x],nrow(transcript_ranges)),exonLen))
}

mm9_reduced <- cbind(gtf_reduced[,c(1:3,5)],1:nrow(gtf_reduced),gtf_reduced[,c(6:7)])
write.table(mm9_reduced,'/media/shaoqizhu/easystore/Annotation/mm9_2015_reduced.bed',quote = F,col.names = F,sep='\t',row.names = F)

mm10_reduced <- gtf_reduced[,c(1:3,5:7)]
mm10_2015_reduced <- cbind(gtf_reduced[,c(1:3,5)],1:nrow(gtf_reduced),gtf_reduced[,c(6:7)])
write.table(mm10_2015_reduced,'/media/shaoqizhu/easystore/Annotation/mm10_2015_reduced.bed',quote = F,col.names = F,sep='\t',row.names = F)
promoter <- reduced
for (i in 1:nrow(reduced)){
  if(reduced$strand[i]=='+'){
    ifelse(reduced$start[i]>1000,
           promoter$start[i] <- reduced$start[i]-1000,
           promoter$start[i] <- 0)
    promoter$end[i] <- reduced$start[i]+1000
  }
  if(reduced$strand[i]=='-'){
    ifelse(reduced$end[i]>1000,
           promoter$start[i] <- reduced$end[i]-1000,
           promoter$start[i] <- 0)
    promoter$end[i] <- reduced$end[i]+1000
  }
}



mm9 <- as.data.frame(cbind(seqnames = as.character(aggregate(gtf$seqnames,list(gtf$gene_name),function(x)x[1])[,2]),
             start = as.character(aggregate(gtf$start,list(gtf$gene_name),min)[,2]),
             end = as.character(aggregate(gtf$end,list(gtf$gene_name),max)[,2]),
             gene_name = sort(unique(gtf$gene_name)),
             strand = as.character(aggregate(gtf$strand,list(gtf$gene_name),function(x)x[1])[,2])))
write.table(mm10,'/media/shaoqizhu/easystore/Annotation/mm10_2015.bed',quote = F,col.names = F,sep='\t',row.names = F)

exon <- gtf[gtf$type=='exon',]
exonLen <- aggregate(exon$width,list(exon$gene_name),sum)
exonLen <- exonLen[order(exonLen$Group.1),]


gencode <- as.data.frame(gtf[gtf$type=='gene'])
gencode <- as.data.frame(gtf[gtf$type=='gene'&gtf$gene_type=='protein_coding'])
gencode <- gencode[,c("seqnames","start","end","gene_name","gene_id","strand")]
write.table(gencode,'~/Documents/annotation/gencode_mm10.bed',quote = F,col.names = F,row.names = F,sep = '\t')


#args = commandArgs(trailingOnly=TRUE)
args1 = '/media/shaoqizhu/easystore/LSD1/mm10_2015.gtf'
# args2 = 'gencode'
# args3 = 'mm9' 
# args4 = '~/Documents/annotation/'

if(length(table(table(gencode$gene_name)))){
  gencode <- gencode[order(gencode$gene_name),]
  gencode.rep <- as.data.frame(table(gencode$gene_name)[table(gencode$gene_name)>1])
  idx=as.numeric(); for(i in 1:nrow(gencode.rep)){idx = c(idx,0:(gencode.rep$Freq[i]-1))}
  rep=gencode$gene_name %in% gencode.rep$Var1
  gencode$gene_name[rep][idx!=0]=paste(gencode$gene_name[rep][idx!=0],idx[idx!=0],sep='-',)
}

if (args2%in%c('gencode','ensembl')){
    geneAnno <- geneAnno[,c("gene_name","seqnames","start","end","strand")]
    if(length(table(table(geneAnno$gene_name)))){
      geneAnno <- geneAnno[order(geneAnno$gene_name),]
      geneAnno.rep <- as.data.frame(table(geneAnno$gene_name)[table(geneAnno$gene_name)>1])
      idx=as.numeric(); for(i in 1:nrow(geneAnno.rep)){idx = c(idx,0:(geneAnno.rep$Freq[i]-1))}
      rep=geneAnno$gene_name %in% geneAnno.rep$Var1
      geneAnno$gene_name[rep][idx!=0]=paste(geneAnno$gene_name[rep][idx!=0],idx[idx!=0],sep='-')
    }
}
refseq <- as.data.frame(gtf[gtf$type=='gene'])
refseq <- refseq[,c("seqnames","start","end","strand","Name")]

if (args2=='refseq'){
    refseq <- refseq[,c("seqnames","start","end","strand","Name")]
    annotation <- as.data.frame(gtf)
    annotation <- annotation[annotation$genome %in% c('chromosome','mitochondrion','genomic') &
                               annotation$type=='region',c('seqnames','width','Name')]
    #write_delim(annotation,paste(args4,'/',args3,'.genome',sep=''),delim = "\t",col_names=F)
    if(length(table(table(refseq$Name)))){
      refseq <- refseq[order(refseq$Name),]
      refseq.rep <- as.data.frame(table(refseq$Name)[table(refseq$Name)>1])
      idx=as.numeric(); for(i in 1:nrow(refseq.rep)){idx = c(idx,0:(refseq.rep$Freq[i]-1))}
      rep=refseq$Name %in% refseq.rep$Var1
      refseq$Name[rep][idx!=0]=paste(refseq$Name[rep][idx!=0],idx[idx!=0],sep='-')
    }
}

geneAnno <- read.table('~/Downloads/fasta/geneAnno.bed')
length(unique(gencode$gene_name))
length(intersect(gencode$gene_name,geneAnno$V5))
gencode <- gencode[order(gencode$gene_name),]
gencode <- gencode[gencode$gene_name %in% intersect(gencode$gene_name,geneAnno$V5),]
write.table(gencode,'~/Downloads/fasta/gencode.bed',quote = F, sep = "\t",row.names = F,col.names = F)

length(unique(refseq$Name))
refseq <- refseq[order(refseq$Name),]
refseq <- refseq[refseq$Name %in% intersect(refseq$Name,geneAnno$V5),]
refseq <- refseq[which(str_split_fixed(refseq$seqnames,'_',2)[,1]=="NC"),]

write.table(refseq,'~/Downloads/fasta/refseq.bed',quote = F, sep = "\t",row.names = F,col.names = F)


gtf <- import(args[1])
geneAnno <- as.data.frame(gtf[gtf$type=='exon'])
gene <- as.data.frame(gtf[gtf$type=='gene'])

exonLen <- aggregate(geneAnno$width,list(geneAnno$gene_name),sum)
exonLen <- exonLen[order(exonLen$Group.1),]
setdiff(gene$gene_name,exonLen$Group.1)

length(unique(gene$gene_name))


mm10 <- read.table('/media/shaoqizhu/easystore/Annotation/mm10_2015.bed')
