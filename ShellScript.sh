#deeptools计算heatmap矩阵
~/.local/bin/computeMatrix reference-point --referencePoint TSS -b 3000 -a 3000 -R /media/shaoqi-zhu/work/ChipSeq/Annotation/TAIR_genes.bed -S /media/shaoqi-zhu/work/ChipSeq/Chip/m5.bw --skipZeros -o
/media/shaoqi-zhu/work/ChipSeq/Chip/m5Matrix.gz
 
#热图
~/.local/bin/plotHeatmap -m matrix_m5.gz -out Heatmap_m5.png
 
#gff转bed文件
cat /media/shaoqi
zhu/work/ChipSeq/Annotation/TAIR10_GFF3_genes_transposons.gff | grep transposable_element_gene | cut -f1,4,5,9 |cut -c4-|cut -f1 -d";"| awk '{print $1, $2, $3, $4}' | sed -e 's/ /\t/g' | sed -e 's/\"//g' > TAIR_TEG.bed
 
#取两个bed文件的交集(-d)和余集(-u)
TAIR_TEA.bed TAIR_TEG.bed | sort | uniq -u > TAIR_TE.bed

#CPU内存状况
gnome-system-monitor
 
#并行计算bigwig文件
for i in {"Col","m5","m6","M51","M52","M61","M62"};do ~/.local/bin/bamCompare -b1 /media/shaoqi-zhu/work/ChipSeq/Chip/${i}-2.bam -b2 /media/shaoqi-zhu/work/ChipSeq/Chip/${i}-1.bam -o /media/shaoqi-zhu/work/ChipSeq/Chip/${i}.bed & done
 
#Cellranger 比对单细胞序列
export PATH=/media/shaoqi-zhu/work/Heart/reference/cellranger-3.0.2:$PATH
cellranger count --id=heart --fastq=/media/shaoqi-zhu/work/Heart/ --transcriptome=/media/shaoqi-zhu/work/Heart/reference/refdata/ --sample=Undetermined,heart1,heart2,heart3,heart4
 
#SRA RNA-seq处理流程
#设置环境变量
export PATH=$PATH:/Users/zsq/samtools
export PATH=$PATH:/Users/zsq/sratoolkit/bin
export PATH=$PATH:/Users/zsq/hisat2
永久添加环境变量：vi  ~/.bash_profile  esc 输入:wq 退出

#从服务器下载文件
scp szhu@tang.phys.gwu.edu:/home/szhu/ncbi/public/sra/SRR3666772_sorted.bam /Users/zsq/ncbi/RP
 
#终止进程
ps -ef
kill -s 9 1827
#批量删除
rm test-{1..20}.txt
#查看nohup.out日志
tail -fn 50 nohup.out

#DogFinder
export PATH=$PATH:~/package/bin
nohup Pre_Process -Q 1 -bam SRR3666751_sorted.bam -ref ~/ref/mm9_loci_annotation.bed &
for i in {SSR3666751,SSR3666755,SSR5558714,SSR5558720}; do nohup Get_DoGs -out ~/rt -bam ${i}_sorted.sorted_PE.bam -suff ${i} -a ~/ref/mm10_loci_annotation.bed & done

#Trimming and chipseq
for i in {"CD4_Tcf1","HSCs_Tle3","LSCs_Tle3","Treg_EED","Treg_EZH2","Treg_HDAC1","p45_Tcf1"}; 

fastqc -o ~/chipseq/trim HSCs_Tle3_R1_val_1.117bp_3prime.fq.gz HSCs_Tle3_R1_val_2.117bp_3prime.fq.gz

trim_galore --hardtrim5 120 --length 30 --cores 4 --paired ${i}_R1.fq.gz ${i}_R2.fq.gz --gzip -o ~/chipseq/trim2 & done
trim_galore --stringency 3 --length 30 --fastqc --paired ${i}_R1.fastq.gz ${i}_R2.fastq.gz --gzip -o ~/chipseq/trim
bowtie2 --mm -p 8 --no-unal --no-mixed --non-deterministic --no-discordant -x ~/ncbi/public/index/bowtie2index/mm10 -1 ${i}_R1.fq.gz -2 ${i}_R2.fq.gz -S ~/chipseq/trim/${i}_trimmed.sam & done
bowtie2 -p 10 -x ~/ncbi/public/index/bowtie2index/mm10 -1 ${i}_R1.fq.gz -2 ${i}_R2.fq.gz | samtools sort -O bam -@ 10 -o -> ${i}.bam

tar -xzvf mm10_library.tar.gz

#序列比对
prefetch
fastq-dump --gzip --split-3 -A SRR366675${i}.sra &
hisat2 -t -mm -p 10 -x ~/ncbi/public/index/mm9/mm9 -1 SRR366675${i}.sra_1.fastq.gz -2 SRR366675${i}.sra_2.fastq.gz -S SRR366675${i}.sam & done
hisat2 -t -mm -p 10 -x ~/ncbi/public/index/hisat2index/grch38/genome -U SRR2033047.sra.fastq.gz -S SRR2033047.sam &

samtools view -h -q 10 -b -f 2 -F 4 -F 8 -F 256 -F 512 -F 2048 ${i}_trimmed.sam > ${i}_trimmed.bam & done
samtools view -@ 10 -b -h -q 10 SRR2033047_grch38.sam > SRR2033047_grch38.bam
samtools sort ${i}_trimmed.bam -o ${i}_sorted.bam & done
samtools index ${i}_sorted.bam & done

#bamCoverage
bamCoverage -b ${i}_sorted_rmdup.bam --binSize 30 --smoothLength 100 --normalizeUsing RPKM -p 10 -o ${i}.bw & done

#Picard deduplication
java -jar ~/package/picard.jar MarkDuplicates I=${i}_sorted.bam O=${i}_sorted_rmdup.bam ASSUME_SORTED=TRUE M=${i}_marked_dup_metrics.txt REMOVE_DUPLICATES=true CREATE_INDEX=true QUIET=true & done

#查看某个染色体
samtools view -b CD4_ATAC_rmdup.bam chrM > CD4_ATAC_rmdup_mitochondria.bam

Find integenic region
gzcat Mus_musculus.GRCm38.97.gtf.gz | grep -v '^#' | awk '$3="gene" {print $1,$4,$5,$14}' | grep -v ENSMUST | sed 's/\"//g' | sed 's/\;//g' | grep -v Gm | wc -l

gzcat gencode.vM22.annotation.gtf.gz | awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4,$5,$12,$14}' | sed 's/\"//g' | sed 's/\;//g' | grep protein_coding | gzip > gencode_gene.bed.gz
bedtools sort -i gencode_gene.bed.gz -g mouse.mm10.genome | bedtools complement -i stdin -g mouse.mm10.genome | gzip > gencode_intergenic.bed.gz

bedtools makewindows -b gencode_intergenic_sorted.bed -w 100 -i srcwinnum | sed 's/\"//g' | awk -F '.' '{print $1"\t"$2}' | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$4}' | gzip > intergenic_windows_order.bed.gz

bedtools sort -i intergenic_windows_order.bed.gz -g mouse.mm10.genome | gzip > intergenic_windows_sorted.bed.gz

bedtools coverage -sorted -g mouse.mm10.genome -a intergenic_windows_sorted.bed.gz -b KCl_sorted.bam | gzip > intergenic_KCl.bed.coverage.gz

gzcat intergenic_KCl.bed.coverage.gz | cut -f4,6 > reads_count.bed.coverage

tar -xvzf Homo_sapiens_NCBI_GRCh38.tar.gz

for i in *than*; do cat $i | awk 'BEGIN{OFS="\t"} {if($2>50000) $2=$2-50000; else $2=0} {$3=$3+50000} {print$0}' > $(cut -d'.' -f1<<<$i)_50kb.bed; done



#行首或行末添加字符
sed 's/^/chr&/g' sed 's/$/&TAIL/g'


gzcat KCl_sorted.bed.gz | grep '\-$' | bedtools coverage -a intergenic_windows_Tnp2.bed.gz -b stdin

samtools sort -n -@ 8 KCl_sorted.bam > KCl_sorted_byname.bam

bedtools bamtobed -i KCl_sorted_byname.bam -bedpe | gzip > KCl_sorted_pe.bed.gz

paste <(gzcat KCl_sorted_pe.bed.gz | cut -f1,2,6) <(gzcat KCL_sorted_r1.bed.gz) | gzip > KCl_sorted_pe_strand_test.bed.gz

#macs2 callpeak
macs2 callpeak -t ${i} -f BAMPE -g mm -n ${i} -B -q 0.1 --SPMR --outdir ~/Downloads/Bigwig/macs2

for i in *pileup.bdg; do bdg2bw $i ~/hdac/bigWig/chromInfo.txt & done

#Homer motif analysis
findMotifsGenome.pl chipseq.bed mm10 ~/Downloads/BigWig/macs2/homer -size 300 -mask

#sort by the 9th column numerically by reverse order
sort -k 9nr chipseq_no_cutrun.bed

cat chipseq_summits.bed | awk 'BEGIN{OFS="\t"} {$2=$2-149}{$3=$3+149} {print $0}' > chipseq.bed


gzcat GCF_000001405.39_GRCh38.p13_genomic.gff.gz | grep -v "^#" | grep NC_ | grep ID=gene | cut -f1 -d ";" | sed 's/ID=gene-//g' | cut -f1,4,5,7,9 > gene.bed


link=$(cat assembly_summary_refseq.txt | grep 'Mus musculus' | grep Chromosome | grep Patch | cut -f20)
file=$(echo $link | rev | cut -d'/' -f1 | rev)
wget $link/$file"_genomic.gtf.gz"


gzcat GCF_000001405.39_GRCh38.p13_genomic.gff.gz | grep 'sequence-region NC' | awk '{print $2"\t"$4}' > grch38.genome



paste -d "-" <(paste -d ":" <(cat geneAnno.bed | cut -f1) <(cat geneAnno.bed | awk '$2=$3-199 {print $2}')) <(cat geneAnno.bed | cut -f3) > polyA_ensembl.bed

paste -d "-" <(paste -d ":" <(cat gencode.bed | grep -w + | cut -f1) <(cat gencode.bed | grep -w + | awk '$7=$3-99 {print $7}')) <(cat gencode.bed | grep -w + | awk '$8=$3+100 {print $8}') > polyA_plus.bed

paste -d "-" <(paste -d ":" <(cat gencode.bed | grep -w - | cut -f1) <(cat gencode.bed | grep -w - | awk '$7=$2-99 {print $7}')) <(cat gencode.bed | grep -w - | awk '$8=$2+100 {print $8}') > polyA_minus.bed

paste -d "-" <(paste -d ":" <(cat gencode.bed | grep -w + | cut -f1) <(cat gencode.bed | grep -w + | awk '$7=$3-199 {print $7}')) <(cat gencode.bed | grep -w + | cut -f3) > polyA_plus.bed

paste -d "-" <(paste -d ":" <(cat gencode.bed | grep -w - | cut -f1) <(cat gencode.bed | grep -w - | cut -f2)) <(cat gencode.bed | grep -w - | awk '$8=$2+199 {print $8}') > polyA_minus.bed

samtools faidx GRCm38.p6.genome.fa -n 200 -i -r polyA_minus.bed > polyA_minus.fa
samtools faidx GRCm38.p6.genome.fa -n 200 -r polyA_plus.bed > polyA_plus.fa

./RNAshapes -f polyA_plus.fa -s -c 5 -t 1 -w 100 -W 100 -O 'D{%s\n}' > polyA_plus_1.txt

python combining_substructure.py -i polyA_plus_1.txt -o polyA_plus_2.txt

python filtering_number_of_ss.py -n 3 -i polyA_plus_2.txt -o polyA_plus_3.txt

python shape_assign_per_nucleotide.py -c 3 -i polyA_plus_3.txt -o polyA_plus_4.txt

python DeepPASTA_polyA_site_prediction_testing.py -testSeq polyA_plus.fa -testSS polyA_plus_4.txt -o polyA_plus.txt

sed -i 's/rc/rc_1/g' polyA_minus.fa

./RNAshapes -f polyA_minus.fa -s -c 5 -t 1 -w 100 -W 100 -O 'D{%s\n}' > polyA_minus_1.txt

python combining_substructure.py -i polyA_minus_1.txt -o polyA_minus_2.txt

python filtering_number_of_ss.py -n 3 -i polyA_minus_2.txt -o polyA_minus_3.txt

python shape_assign_per_nucleotide.py -c 3 -i polyA_minus_3.txt -o polyA_minus_4.txt

python DeepPASTA_polyA_site_prediction_testing.py -testSeq polyA_minus.fa -testSS polyA_minus_4.txt -o polyA_minus.txt


cat na-WT_D0.bed | awk '{if($5>5 || $5<-2.5) print$0}'>na-WT_D0_cut.bed
bedtools intersect -a <(cat na-WT_D0.bed | awk 'BEGIN{OFS="\t";} $5>1 {print$0}') -b <(cat dKO_D0-WT_D0.bed | awk 'BEGIN{OFS="\t";} $5>1 {print$0}') > WT_D0-na_dKO_D0.bed
bedtools intersect -a <(cat WT_D0-WT_D1.bed | awk 'BEGIN{OFS="\t";} $5>1 {print$0}') -b <(cat dKO_D1-WT_D1.bed | awk 'BEGIN{OFS="\t";} $5>1 {print$0}') > WT_D1-WT_D0_dKO_D1.bed


#linux library update for r package installation
curl: sudo apt-get install curl
libssl-dev: sudo apt-get install libssl-dev
libcurl: sudo apt-get install libcurl4-openssl-dev
xml2: sudo apt-get install libxml2-dev


