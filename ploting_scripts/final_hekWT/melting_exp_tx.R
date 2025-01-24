
############################################################################
# Load libraries
library(ggplot2)
library(dplyr)
library(reshape2)
library(pheatmap)
library(ComplexHeatmap)
library("dendextend") 
library(ggsci)
library(eulerr)
library(tidyr)
library(biomaRt)
mypal <- pal_npg("nrc", alpha = 0.7)(9)

rna_hek<-read.table('hek293_rnaseq_S1_rsem.isoforms.results',sep = '\t',header = T)
head(rna_hek)
rna_hek$transcript_id_clean <- sub("\\..*", "", rna_hek$transcript_id)

convertible_df<-rna_hek

head(convertible_df)
cor(convertible_df$expected_count,convertible_df$length)

rna<-read.table('hg19_rna_seq.bed.gz',sep = '\t',header = T)
long_gene_raw<-rna[(rna$end-rna$start)>300000,]
colnames(long_gene_raw)<-c('chr','start','end','transcript','type','signal','gene','info')
long_gene_raw$transcript_id_clean<-sub("\\..*", "", long_gene_raw$transcript)
long_gene_hek_exp <- merge(long_gene_raw, convertible_df, by = "transcript_id_clean")
colnames(long_gene_hek_exp)[14]
#long_gene_12878_exp<-long_gene_12878_exp[long_gene_12878_exp$TPM_ci_upper_bound>0.05,]
long_gene<-long_gene_hek_exp[,c(2,3,4,5,6,14,8,9)]

colnames(long_gene)<-c('chr','start','end','transcript','type','signal','gene','info')
long_gene$signal<-log10(long_gene$signal+1)
melt_sum<-read.table('tx_melt_score.xls')
colnames(melt_sum)<-c('transcript','ins_mean_clust_1','ins_mean_clust_2','d','pval','cellNum_clust_1','cellNum_clust_2')

melt_sum$clust_2_clust_1_ratio<-melt_sum$cellNum_clust_1/(melt_sum$cellNum_clust_1+melt_sum$cellNum_clust_2)

melt_sum_exp<-merge(melt_sum,long_gene,by='transcript')
melt_sum_exp$exp_norm_by_length<-melt_sum_exp$signal/(melt_sum_exp$end-melt_sum_exp$start)
melt_sum_exp<-melt_sum_exp[melt_sum_exp$d<1,]

ggplot(melt_sum_exp, aes(x = melt_sum_exp$signal)) +  geom_density()
ggplot(melt_sum_exp, aes(x = log10(melt_sum_exp$signal))) +  geom_density()

#filtering a/b
library(data.table)
eigen<-fread('snHiChew_HEK293Twt.mcool.res_500000.eig.cis.vecs.tsv')
eigen_b<-eigen[eigen$E1<0,]
colnames(eigen_b)[1]<-'chr'
eigen_b_df<-as.data.frame(eigen_b)

is_within_range <- function(chr, start, end, table2) {
  chr_table2 <- table2 %>% filter(chr == !!chr)
  any(apply(chr_table2, 1, function(range_row) {
    (end >= as.numeric(range_row['start']) && end <= as.numeric(range_row['end'])) ||
      (start >= as.numeric(range_row['start']) && start <= as.numeric(range_row['end']))  }))}

thres<-2
# Apply the filter
filtered_melt_sum_exp <- melt_sum_exp %>%  rowwise() %>%  filter(is_within_range(chr, start, end, eigen_b))
b_comp_high<-length(rownames(filtered_melt_sum_exp[filtered_melt_sum_exp$signal>thres,]))
b_comp_low<-length(rownames(filtered_melt_sum_exp[filtered_melt_sum_exp$signal<=thres,]))

melt_sum_exp_a<-melt_sum_exp[!melt_sum_exp$transcript %in% filtered_melt_sum_exp$transcript,]
a_comp_high<-length(rownames(melt_sum_exp_a[melt_sum_exp_a$signal>thres,]))
a_comp_low<-length(rownames(melt_sum_exp_a[melt_sum_exp_a$signal<=thres,]))

# Create a sample data frame
data <- data.frame( category = c("b_comp_high", "b_comp_low"),
  count = c(b_comp_high,b_comp_low))

# Calculate the percentage for each category
data$percentage <- data$count / sum(data$count) * 100

p0<-ggplot(data, aes(x = "", y = percentage, fill = category)) +
  geom_bar(width = 1, stat = "identity",alpha=0.7) +
  coord_polar("y", start = 0) +
  theme_void() + # Remove background, grid, and numeric labels
  labs(title = paste0('  b_comp_high:',b_comp_high,' b_comp_low:',b_comp_low)) +
  theme(legend.title = element_blank())+scale_fill_manual(values=c(mypal[1],mypal[3]))

p0
filtered_melt_sum_exp$exp_level<-'transcript_low_exp'
filtered_melt_sum_exp[filtered_melt_sum_exp$signal>thres,length(colnames(filtered_melt_sum_exp))] <-'transcript_high_exp'

melt_sum_exp_high<-filtered_melt_sum_exp[filtered_melt_sum_exp$signal>thres,]
melt_sum_exp_low<-filtered_melt_sum_exp[filtered_melt_sum_exp$signal<thres,]

#write.table(melt_sum_exp_high, file = 'melt_sum_exp_high.xls',sep = "\t",col.names=T,row.names=FALSE,quote=FALSE)
round(cor(log2(melt_sum_exp_high$end-melt_sum_exp_high$start),melt_sum_exp_high$d),3)
round(cor(log2(melt_sum_exp_high$end-melt_sum_exp_high$start),melt_sum_exp_high$signal),3)
round(cor(log2(melt_sum_exp_high$d),melt_sum_exp_high$signal),3)

melt_sum_exp_high$d
thres
#high
pcc<-round(cor(log2(melt_sum_exp_high$exp_norm_by_length),melt_sum_exp_high$d),3)
geneNum<-length(melt_sum_exp_high$exp_norm_by_length)
p1<-ggplot(melt_sum_exp_high, aes(x = melt_sum_exp_high$d,y=log2(melt_sum_exp_high$exp_norm_by_length))) +  
  geom_point(alpha=0.5,size=3,stroke=0,color='red3')+theme_classic()+geom_smooth(method='lm',color='grey20',linetype='dashed')+
  labs(title=paste0('PCC:',pcc,', transcript Num:',geneNum),subtitle = paste0('high exp long transcript (log10 RSEM Expected Count > ',thres,')'),
       x='KS test maximum distance',y='log10 (RSEM Expected Count) length scaled transcript exp')
p1

wilcox_p<-wilcox.test(melt_sum_exp_low$d,melt_sum_exp_high$d)

p2<-ggplot(filtered_melt_sum_exp, aes(x = d,group=exp_level,color=exp_level)) +  
  geom_density(lwd=1)+theme_classic()+scale_color_manual(values=c(mypal[1],mypal[3]))+
  labs(title=paste0(subtitle = paste0('wilcox p: ',round(as.numeric(wilcox_p[3]),3),
  '\nhigh exp long transcript (log10 RSEM Expected Count > ',thres,')')),x='KS test maximum distance',y='density')
p2
#ggplot(filtered_melt_sum_exp, aes(y = d,x=exp_level,color=exp_level)) + geom_violin()
#low
pcc<-round(cor(log2(melt_sum_exp_low$exp_norm_by_length+0.00001),melt_sum_exp_low$d),3)
geneNum<-length(melt_sum_exp_low$exp_norm_by_length)
p3<-ggplot(melt_sum_exp_low, aes(x = melt_sum_exp_low$d,y=log2(melt_sum_exp_low$exp_norm_by_length))) +  
  geom_point(alpha=0.5,size=3,stroke=0,color=mypal[3])+theme_classic()+geom_smooth(method='lm',color='grey20',linetype='dashed')+
  labs(title=paste0('PCC:',pcc,', transcript Num:',geneNum),subtitle = paste0('low exp long transcript (log10 RSEM Expected Count < ',thres,')'),
       x='KS test maximum distance',y='log10 (RSEM Expected Count) length scaled transcript exp')
p3




library(gridExtra)

pdf('melting_b.pdf',width = 11,height=14)
grid.arrange(p0,p2,p1,p3,nrow=2)
dev.off()

wt_kd_RNA_meltB_high$Software_Fold_Change
####################################################################################################
wt_kd_RNA<-read.table('groupA-vs-groupC.gene.DEseq2.xls',sep = '\t',header = T)
#wt_kd_RNA<-read.table('groupA-vs-groupC.gene.DEseq2.Filter.xls',sep = '\t',header = T)

wt_kd_RNA_meltB_high<-merge(melt_sum_exp_high,wt_kd_RNA,by.x='gene',by.y='GeneID')
wt_kd_RNA_meltB_high$group<-'sig_fc1.5'
wt_kd_RNA_meltB_high[wt_kd_RNA_meltB_high$Corrected_Significant_Difference_Value>0.05|abs(wt_kd_RNA_meltB_high$Software_Fold_Change)<0.58,61]<-'non_sig'



p7<-ggplot(wt_kd_RNA_meltB_high, aes(x=wt_kd_RNA_meltB_high$Software_Fold_Change,y= -log10(wt_kd_RNA_meltB_high$Corrected_Significant_Difference_Value),color=group)) + geom_point(alpha=1,size=1.5)+ 
  scale_color_manual(values=c("grey",'#673147'))+geom_hline(yintercept=c(1.30103),linetype="dotted")+geom_vline(xintercept=c(0.58,-0.58),linetype="dotted")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold")) + xlim(-2,2)  +
  labs(title=paste0('melting b long high ',length(rownames(wt_kd_RNA_meltB_high)),' genes\nX: log2(kd / wt), Y: -log10(padj)
                    \n#sig_fc1.5: ',length(unique(wt_kd_RNA_meltB_high[wt_kd_RNA_meltB_high$group=='sig_fc1.5',1]))))
p7

#melt_sum_exp_low
wt_kd_RNA_meltB_low<-merge(melt_sum_exp_low,wt_kd_RNA,by.x='gene',by.y='GeneID')
wt_kd_RNA_meltB_low$group<-'sig_fc1.5'
wt_kd_RNA_meltB_low[wt_kd_RNA_meltB_low$Corrected_Significant_Difference_Value>0.05|abs(wt_kd_RNA_meltB_low$Software_Fold_Change)<0.58,61]<-'non_sig'
length(unique(wt_kd_RNA_meltB_low[wt_kd_RNA_meltB_low$group=='sig_fc1.5',1]))

p8<-ggplot(wt_kd_RNA_meltB_low, aes(x=wt_kd_RNA_meltB_low$Software_Fold_Change,y= -log10(wt_kd_RNA_meltB_low$Corrected_Significant_Difference_Value),color=group)) + geom_point(alpha=1,size=1.5)+ 
  scale_color_manual(values=c("grey",'#673147'))+geom_hline(yintercept=c(1.30103),linetype="dotted")+geom_vline(xintercept=c(0.58,-0.58),linetype="dotted")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold")) + xlim(-2,2)  +
  labs(title=paste0('melting b long low ',length(rownames(wt_kd_RNA_meltB_low)),' genes\nX: log2(kd / wt), Y: -log10(padj)
                    \n#sig_fc1.5: ',length(unique(wt_kd_RNA_meltB_low[wt_kd_RNA_meltB_low$group=='sig_fc1.5',1]))))
p8
dev.off()
pdf('melting_b_long_gene_bulkRNA.pdf',width = 12,height=5)
grid.arrange(p7,p8,nrow=1)
dev.off()

