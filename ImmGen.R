library(openxlsx)
cell_type <- c(rep('splenic_B',3),rep('slpenic_DC',8),rep('BM_granulocytes',2),rep('splenic_NK',5),
               rep('splenic_naive_CD4',2),rep('splenic_naive_CD8',2),rep('CD8_effectors',2),
               rep('splenic_gd_T',2),rep('Treg_cells',2),rep('monocytes_in_PBLs',4))
cell_idx <- c(rep(1,3),rep(2,8),rep(3,2),rep(4,5),rep(5,2),rep(6,2),rep(7,2),rep(8,2),rep(9,2),rep(10,4))
cell <- as.data.frame(cbind(cell_type,cell_idx))

immune_cell <- read.xlsx('/media/shaoqizhu/easystore/CD8_HP_stimulated/CD8_HP_naive/GSE109125_selected subsets.xlsx')
#which(immune_cell$gene_symbol==43526)
#immune_cell <- immune_cell[-c(39039,39040),]
#immune_cell <- immune_cell[-c(39041,39044),]
immune_cell[39041,1] <- '43525-1'
immune_cell[39044,1] <- '43526-1'
rownames(immune_cell) <- immune_cell$gene_symbol
immune_cell <- immune_cell[,-1]

immune_cell <- as.data.frame(t(t(immune_cell)/apply(immune_cell, 2, sum)*mean(apply(immune_cell,2,sum))))

for(i in 1:10){
  thres = 1
  cell.1 = which(cell$cell_idx==i)
  cell.2 = which(cell$cell_idx!=i)
  data.1 = apply(X = immune_cell[, cell.1, drop = F], MARGIN = 1, 
                 FUN = function(x) log(x = mean(x)))
  data.2 = apply(X = immune_cell[, cell.2, drop = F], MARGIN = 1, 
                 FUN = function(x) log(x = mean(x)))
  avg_logFC = data.1-data.2
  genes = names(avg_logFC)[abs(avg_logFC)>thres]
  group.info <- data.frame(row.names = c(cell.1, cell.2))
  group.info[cell.1, "group"] <- "Group1"
  group.info[cell.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(x = group.info[, "group"])
  p_val <- pbsapply(
    X = genes,FUN = function(x) {
      return(t.test(t(log(immune_cell[x,])) ~ group.info[, "group"])$p.value)})
  p_val_adj <- p.adjust(p = p_val,method = 'BH',n = length(genes))
  
  marker = as.data.frame(cbind(avg_logFC = signif(avg_logFC[genes],7),
                            p_val = signif(unlist(p_val),7),
                            p_val_adj = signif(p_val_adj,7)))
  rownames(marker) <- genes
  DEG <- marker[marker$p_val_adj<0.01 & abs(marker$avg_logFC)>1,]
  write.csv(DEG[order(DEG$avg_logFC,decreasing = T),],paste('~/Documents/',unique(cell$cell_type[cell.1]),'.csv',sep=''))
}
cell$cell_type[cell.1]
unique(cell$cell_type[cell.1])
immune_cell['Cd5',cell.1]

immune_cell <- read.csv('/media/shaoqizhu/easystore/CD8-HP/CD8_HP_naive/GSE109125_Selected_Subsets_Gene_Count_31.txt',sep = '\t')
rownames(immune_cell) <- immune_cell$Gene_Symbol; immune_cell <- immune_cell[,-1]
cell_type <- as.data.frame(cbind(names(immune_cell),c(rep('Bcell',2),rep('mDC',6),rep('pDC',2),rep('Granulocyte',2),rep('NK',5),rep('Tcell',10),rep('Monocyte',4))))

Tcell <- cell_type[cell_type$V2=='Tcell',]
Tcell$V2 <- c(rep('Teff',2),rep('CD4',2),rep('CD8',2),rep('Tgd',2),rep('Treg',2))
immune_Tcell <- immune_cell[,Tcell$V1]

diff <- numeric(); N = length(unique(cell_type$V2))
for (i in 1:(N-1)){
  for (j in (i+1):N){
    condition1 = which(cell_type$V2==unique(cell_type$V2)[i])
    condition2 = which(cell_type$V2==unique(cell_type$V2)[j])
    group <- factor(c(rep(i,length(condition1)),rep(j,length(condition2))))
    y = DGEList(counts = immune_cell[,c(condition1,condition2)],group = group)
    print(names(immune_cell[,c(condition1,condition2)]))
    y <- calcNormFactors(y)
    y <- estimateCommonDisp(y)
    fit <- glmQLFit(y)
    qlf <- glmQLFTest(fit)
    fdr <- p.adjust(p = qlf$table$PValue,method = "BH",n = length(qlf$table$PValue))
    diff_idx <- which(fdr < 0.05 & abs(qlf$table$logFC) > log2(2))
    diff <- c(diff,diff_idx)
    change <- cbind(Gene_Symbol=rownames(immune_cell)[diff_idx],logFC=qlf$table$logFC[diff_idx],pval=qlf$table$PValue[diff_idx],fdr=fdr[diff_idx])
    write.table(change,paste('/media/shaoqizhu/easystore/CD8-HP/CD8_HP_naive/',unique(cell_type$V2)[i],'-',unique(cell_type$V2)[j],'.bed',sep=''),
               quote = F,sep = '\t',row.names = F,col.names = T)
  }
}
unique(diff)


Tcell <- cell_type[cell_type$V2=='Tcell',]
Tcell$V2 <- c(rep('Teff',2),rep('CD4',2),rep('CD8',2),rep('Tgd',2),rep('Treg',2))
Tcell$V3 <- c(rep('Teff',2),rep('Tna',8))
immune_Tcell <- immune_cell[,names(immune_cell)%in%Tcell$V1]

diff <- numeric(); N = length(unique(Tcell$V3))
for (i in 1:(N-1)){
  for (j in (i+1):N){
    condition1 = which(Tcell$V3==unique(Tcell$V3)[i])
    condition2 = which(Tcell$V3==unique(Tcell$V3)[j])
    group <- factor(c(rep(i,length(condition1)),rep(j,length(condition2))))
    y = DGEList(counts = immune_Tcell[,c(condition1,condition2)],group = group)
    print(names(immune_Tcell[,c(condition1,condition2)]))
    y <- calcNormFactors(y)
    y <- estimateCommonDisp(y)
    fit <- glmQLFit(y)
    qlf <- glmQLFTest(fit)
    fdr <- p.adjust(p = qlf$table$PValue,method = "BH",n = length(qlf$table$PValue))
    diff_idx <- which(fdr < 0.05 & abs(qlf$table$logFC) > log2(2))
    diff <- c(diff,diff_idx)
    change <- cbind(Gene_Symbol=rownames(immune_Tcell)[diff_idx],logFC=qlf$table$logFC[diff_idx],pval=qlf$table$PValue[diff_idx],fdr=fdr[diff_idx])
    write.table(change,paste('/media/shaoqizhu/easystore/CD8-HP/CD8_HP_naive/Tcell_DEG/',unique(Tcell$V3)[i],'-',unique(Tcell$V3)[j],'.bed',sep=''),
                quote = F,sep = '\t',row.names = F,col.names = T)
  }
}
unique(diff)

CD8HP_naive <- read.table('/media/shaoqizhu/easystore/CD8-HP/CD8_HP_naive/CD8_HP_DKO_vs_WT_naive_read_count.txt',header = T,row.names = 1)
group <- factor(c(rep(1,3),rep(2,3)))
y = DGEList(counts = CD8HP_naive,group = group)
print(names(CD8HP_naive))
y <- calcNormFactors(y,method='TMM')
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
counts_tmm <- t(t(y$pseudo.counts)*y$samples$norm.factors)
write.table(counts_tmm,'/media/shaoqizhu/easystore/CD8-HP/CD8_HP_naive/CD8_HP_DKO_vs_WT_naive_read_count_TMM.txt',quote = F,sep = '\t')

CD8HP_atac_naive <- read.csv('~/Dropbox/Shared_Files/ImmGen_atac_mm10/read_count_mm9_all_ORCs.bed',sep = '\t')

group <- factor(c(1,2,2,1,1))
y = DGEList(counts = CD8HP_atac_naive[,c(5:9)],group = group)
print(names(CD8HP_atac_naive))
y <- calcNormFactors(y,method='TMM')
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
counts_tmm <- t(t(y$pseudo.counts)*y$samples$norm.factors)
write.table(cbind(CD8HP_atac_naive[,1:4],counts_tmm),
            '~/Dropbox/Shared_Files/ImmGen_atac_mm10/read_count_mm9_all_ORCs_TMM.bed',quote = F,sep = '\t')

coldata <- as.data.frame(c(rep(1,3),rep(2,3)))
colnames(coldata) = 'condition'
rownames(coldata) = names(CD8HP_naive)
coldata$condition <- factor(coldata$condition)
dds <- DESeqDataSetFromMatrix(countData = round(CD8HP_naive),
                              colData = coldata,
                              design = ~ condition)
dds <- DESeq(dds,fitType = 'local')
res <- results(dds)
dds


Imm_RNA <- read.csv('/media/shaoqizhu/easystore/CD8-HP/CD8_HP_naive/GSE109125_Selected_Subsets_Gene_Count_31.txt',sep='\t')
Imm_ATAC <- read.csv('/media/shaoqizhu/easystore/ImmuneCell/ImmGenATAC18_AllOCRsInfo.csv')
write.csv(sort(names(Imm_ATAC)),'/media/shaoqizhu/easystore/ImmuneCell/Imm_ATAC_names.csv',row.names = F,col.names = F)
names(Imm_RNA)

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

for (i in unique(condition_ATAC)){
  diff_atac_0 <- diff_atac[,grep(i,colnames(diff_atac))]
  idx <- grep(paste('/',i,sep = ''),colnames(diff_atac_0))
  diff_atac_0[,idx] <- -diff_atac_0[,idx]
min_change <- apply(diff_atac_0,1,min)
print(c(i,length(which(min_change>9))))
}
