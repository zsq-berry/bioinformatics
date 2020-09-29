file_path = '/media/shaoqizhu/easystore/ImmuneCell/'
file_name = 'ImmGenATAC18_AllOCRsInfo.csv'
threshold = 20
Imm_ATAC <- read.csv(paste(file_path,file_name,sep = ''))

name_ATAC <- c("T.8.Nve.Sp","B.Sp","DC.4..Sp","DC.8..Sp","DC.pDC.Sp","GN.BM","NK.27.11b..Sp","NK.27.11b..Sp.1",
               "NK.27.11b..Sp.2","T8.TE.LCMV.d7.Sp","T.4.Nve.Sp","Tgd.Sp","Treg.4.25hi.Sp","Mo.6C.II..Bl","Mo.6C.II..Bl.1")
condition_ATAC <- c("T.8.Nve","B","mDC","mDC","pDC","GN","NK","NK","NK","T8.TE","T.4.Nve","Tgd","Treg","Mo","Mo")

Imm_ATAC_sub <- Imm_ATAC[,name_ATAC]
Imm_ATAC_sub <- as.matrix(Imm_ATAC_sub)
N=length(unique(condition_ATAC)); diff_atac <- as.numeric(); diff_name <- as.character()
for (i in 1:(N-1)){
  for (j in (i+1):N){
    condition1 = which(condition_ATAC==unique(condition_ATAC)[i])
    condition2 = which(condition_ATAC==unique(condition_ATAC)[j])
    print(c(condition1,condition2))
    diff_atac_0 <- unlist(lapply(1:nrow(Imm_ATAC_sub),function(x) mean(Imm_ATAC_sub[x,condition2])-mean(Imm_ATAC_sub[x,condition1])))
    diff_atac <- cbind(diff_atac,diff_atac_0)
    diff_name_0 <- paste(unique(condition_ATAC)[j],'/',unique(condition_ATAC)[i],sep = '')
    diff_name <- c(diff_name,diff_name_0)
  }
}
colnames(diff_atac) <- diff_name

for (i in unique(condition_ATAC)){
  diff_atac_0 <- diff_atac[,grep(i,colnames(diff_atac))]
  idx <- grep(paste('/',i,sep = ''),colnames(diff_atac_0))
  diff_atac_0[,idx] <- -diff_atac_0[,idx]
  min_change <- apply(diff_atac_0,1,min)
  summit <- cbind(index=which(min_change>threshold),min_change=min_change[which(min_change>threshold)],chr=Imm_ATAC$chrom[which(min_change>threshold)],summit=Imm_ATAC$Summit[which(min_change>threshold)],Imm_ATAC_sub[which(min_change>threshold),])
  write.csv(summit,paste(file_path,i,'.csv',sep = ''),row.names = F)
}

