library(ggplot2)
library(plyr)
library(reshape2)
library(zoo)
library(pheatmap)
#library(HiCcompare)
library(Matrix)
library(data.table)
library(scico)
################################################################################################ a/b compartment
hic<-fread("hek293_ngsHiC_700M.mcool.res_250000.eig.cis.vecs.tsv")
hic<-hic[-which(is.na(hic$E1)),]
hic<-hic[hic$chrom=='7',]
hic$locus_id<-rownames(hic)
ggplot(hic,aes(x=rownames(hic),y=hic$E1)) + geom_bar(stat="identity")
hic_line<-rbind.data.frame(cbind.data.frame(locus_id=hic$locus_id,eigen=hic$E1,sample='hic'))
ggplot(hic_line,aes(x=locus_id,y=eigen,group=sample,color=sample)) + geom_line(stat="identity") +
  scale_color_manual(values = c('black')) +theme_bw()+ geom_hline(yintercept = 0,linetype="dashed")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

hic<-fread("GSM3840083_rs.mcool.res_500000.eig.cis.vecs.tsv")
hic<-hic[-which(is.na(hic$E1)),]
sca<-fread("snHiChew_testis_clustLeft.mcool.res_500000.eig.cis.vecs.tsv")
sca<-sca[-which(is.na(sca$E1)),]

sca_hic<-merge(hic,sca,by=c('chrom','start'))
max(sca_hic$GC.x)
min(sca_hic$GC.x)
eigen_diff_gc<-cbind.data.frame(eigen_diff=abs(sca_hic$E1.x-sca_hic$E1.y),GC_ratio=sca_hic$GC.y)
pdf('gc_eigen_diff_curve.pdf')
ggplot(eigen_diff_gc,aes(x=GC_ratio,y=eigen_diff)) + geom_smooth(method = "loess") +
  theme_bw() +scale_color_manual(values=c("purple"))+
  scico::scale_fill_scico(palette = "davos",direction=-1)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+ylim(0,0.5)
 #method = "loess",
dev.off()

sca_hic<-merge(hic,sca,by=c('chrom','start'))
#sca_hic<-sca_hic[sca_hic$GC.x>0.43&sca_hic$GC.x<0.55,]

#sca_hic<-sca_hic[abs(sca_hic$weight.x-sca_hic$weight.y)<0.005,]
cor(sca_hic$E1.x,sca_hic$E1.y)
correlation<-round(cor(sca_hic$E1.x,sca_hic$E1.y),2)

#sca_hic
hic<-hic[hic$chrom=='chr7',]
ggplot(hic,aes(x=rownames(hic),y=hic$E1)) + geom_bar(stat="identity")
sca<-sca[sca$chrom=='chr7',]
ggplot(sca,aes(x=rownames(sca),y=sca$E1)) + geom_bar(stat="identity")

ggplot(sca_hic,aes(x=sca_hic$E1.x,y=sca_hic$E1.y)) + geom_point()
pdf('eigen_snHiChew_testis_clustLeft_GSM3840083_rs.pdf',width=8,height=5)
ggplot(sca_hic, aes(x=sca_hic$E1.x,y=sca_hic$E1.y)) + geom_hex(bins = 80)+ geom_abline(intercept = 0,alpha=0.2, slope = 1,linetype = "dashed") +#geom_smooth(method=lm)+ 
  scico::scale_fill_scico(palette = "davos",direction=-1)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold")) + xlim(-2,2) + ylim(-2,2) +
  labs(title=paste0('number of bins: ',length(rownames(sca_hic))
                    ,"\neigenvector, x:GSM3840083_rs, y:snHiChew_testis_clustLeft, \nmm10, bin 500k, cor ",correlation))
dev.off()
################################################################################################ TAD insulation
hicIns<-fread("GSM3840083_rs.mcool.res_50000.window_100k.ins")
#hicIns<-fread("ngs_hic_2_hek293_hg38.mcool.res_25k.window_100k.ins")
#hicIns<-hicIns[hicIns$chrom=='7',]
hicIns<-hicIns[-which(is.na(hicIns$log2_insulation_score_100000)),]

scaIns<-fread("snHiChew_testis_clustLeft.mcool.res_50000.window_100k.ins")
#scaIns<-scaIns[scaIns$chrom=='7',]
scaIns<-scaIns[-which(is.na(scaIns$log2_insulation_score_100000)),]
#scaIns[scaIns$log2_insulation_score_100000>0,5]<-scaIns[scaIns$log2_insulation_score_100000>0,5]/1.7
sca_hic_Ins<-merge(hicIns,scaIns,by=c('chrom','start'))
head(sca_hic_Ins)
colnames(sca_hic_Ins)[5] <- 'GSM3840083_rs_log2_ins_50k'
colnames(sca_hic_Ins)[10] <- 'snHiChew_testis_clustLeft_log2_ins_50k'

cor(x=sca_hic_Ins$GSM3840083_rs_log2_ins_50k,y=sca_hic_Ins$snHiChew_testis_clustLeft_log2_ins_50k)

#ggplot(sca_hic_Ins,aes(x=rownames(sca_hic_Ins),y=sca_hic_Ins$log2_insulation_score_100000.x)) + geom_bar(stat="identity")
#ggplot(sca_hic_Ins,aes(x=rownames(sca_hic_Ins),y=sca_hic_Ins$log2_insulation_score_100000.y)) + geom_bar(stat="identity")
#ggplot(sca_hic_Ins,aes(x=sca_hic_Ins$log2_insulation_score_100000.x,y=sca_hic_Ins$log2_insulation_score_100000.y)) + geom_point(alpha=0.04) + xlim(-2,2) + ylim(-2,2)

pdf('ins_snHiChew_1200cell_ngsHiC_50k.pdf',width=8,height=5)
ggplot(sca_hic_Ins, aes(x=GSM3840083_rs_log2_ins_50k,y=snHiChew_testis_clustLeft_log2_ins_50k)) + geom_hex(bins = 200)+ geom_abline(intercept = 0,alpha=0.2, slope = 1,linetype = "dashed") + #geom_smooth(method=lm)+
  scico::scale_fill_scico(palette = "davos",direction=-1)+ xlim(-2,2) +ylim(-2,2)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))  +
  labs(title="insulation score, x:HEK_ngsHiC, y:snHiChew_1200cell, hg19, bin 50k, w 100k, cor 0.86")
dev.off()
