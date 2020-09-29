library(HiTC)
library(HiCDataHumanIMR90)
library(BSgenome.Hsapiens.UCSC.hg18)
data(Dixon2012_IMR90)
show(hic_imr90_40)
class(intdata(hic_imr90_40$chr1chr1))
object.size(hic_imr90_40)
seqlevels(hic_imr90_40)
detail(hic_imr90_40$chr6chr6)
head(summary(hic_imr90_40))
sset <- reduce(hic_imr90_40, chr=c("chr5"))
imr90_500 <- HTClist(mclapply(sset, binningC, binsize=50000, bin.adjust=FALSE, method="sum", step=1))
mapC(imr90_500)
resFrag <- getRestrictionFragmentsPerChromosome(resSite="AAGCTT", overhangs5=1,
                                                chromosomes="chr5",
                                                genomePack="BSgenome.Hsapiens.UCSC.hg18")
resFrag
