library(HiCExperiment)
library(rtracklayer)
library(HiContacts)
library(HiContactsData)
library(GenomicInteractions)
library(InteractionSet)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)
library(gridExtra)
library(gtable)
library(ggplot2)
library(BiocParallel)
library(tibble)
library(hicrep)
library(HiCool)

library(ggplot2)
library(plyr)
library(reshape2)
library(zoo)
library(pheatmap)
#library(HiCcompare)
library(Matrix)
library(data.table)
library(scico)
library(gridExtra)
library(patchwork)
library(purrr)

eigen_bin<-'500000'
ins_bin<-'50000'
ins_start<-39000000
ins_end<-50000000

#sample_group<-c('hichew_combine_326M','hichew_SIP_5k','hichew_SIP_500','hichew_SIP_50','hichew_PCR_5k','hichew_PCR_500','hichew_PCR_50')
#sample_group<-c('hichew_SIP_5k','hichew_SIP_500','hichew_SIP_50','hichew_PCR_5k','hichew_PCR_500','hichew_PCR_50')
sample_group<-c('snHiChew_testis_clustLeft')

#sample_list[1]
sample_list<-as.list(sample_group)
#plot_list_eigen <- list()
#plot_list_ins <- list()
#for(i in 1:length(sample_group)){

process_elements <- function(x){
################################################################################################ a/b compartment
sample_1<-x
#sample_1<-'guomei_0428_300M'
sample_2<-'GSM3840083_rs'

eigen_1<-fread(paste0(sample_1,'.mcool.res_',eigen_bin,'.eig.cis.vecs.tsv'))
eigen_2<-fread(paste0(sample_2,'.mcool.res_',eigen_bin,'.eig.cis.vecs.tsv'))

eigen_1<-eigen_1[-which(is.na(eigen_1$E1)),]
eigen_2<-eigen_2[-which(is.na(eigen_2$E1)),]

eigen_merge<-merge(eigen_1,eigen_2,by=c('chrom','start'))
eigen_merge<-eigen_merge[eigen_merge$GC.y>0.4,]
cor(eigen_merge$E1.x,eigen_merge$E1.y,method='spearman')

p_cor<-round(cor(eigen_merge$E1.x,eigen_merge$E1.y,method='pearson'),3)
s_cor<-round(cor(eigen_merge$E1.x,eigen_merge$E1.y,method='spearman'),3)

p1<-ggplot(eigen_merge, aes(x=eigen_merge$E1.x,y=eigen_merge$E1.y)) + geom_hex(bins = 30)+ geom_abline(intercept = 0,alpha=0.2, slope = 1,linetype = "dashed") +#geom_smooth(method=lm)+ 
  scico::scale_fill_scico(palette = "lajolla",direction=-1)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold")) + xlim(-2,2) + ylim(-2,2) +
  theme(legend.direction = "horizontal",legend.position = "top", legend.box = "top")+ labs(x = NULL, y = NULL)+
  labs(title=paste0("eigen    x: ",sample_1,", \ny :",sample_2,", mm10, bin, \npearson: ",p_cor," spearman:",s_cor)) 

pdf(paste0('testis_rou_ab_eigen_cor_20241126.pdf'),width=5,height=5)
p1
dev.off()

eigen_1<-fread(paste0(sample_1,'.mcool.res_',eigen_bin,'.eig.cis.vecs.tsv'))
eigen_2<-fread(paste0(sample_2,'.mcool.res_',eigen_bin,'.eig.cis.vecs.tsv'))
eigen_1<-eigen_1[eigen_1$chrom=='chr7',]
eigen_1$locus_id<-rownames(eigen_1)
eigen_2<-eigen_2[eigen_2$chrom=='chr7',]
eigen_2$locus_id<-rownames(eigen_2)
eigen_1$sample<-sample_1
eigen_2$sample<-sample_2
eigen_line<-rbind.data.frame(eigen_1,eigen_2)

p3<-ggplot(eigen_line,aes(x=locus_id,y=E1,group=sample,color=sample)) + geom_line(stat="identity") +
  scale_color_manual(values = c('#E64B35B2','#4DBBD5B2')) +theme_bw()+ geom_hline(yintercept = 0,linetype="dashed")+
  theme(axis.text.x=element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
  theme(legend.position="top")+ labs(x = NULL, y = NULL)+ theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
p3
#dev.off()

#pdf(paste0('eigen_',sample_1,'_',sample_2,'_',eigen_bin,'.pdf'),width=5,height=10)
#cowplot::plot_grid(p1,p3,nrow=2)


################################################################################################ TAD insulation
ins_1<-fread(paste0(sample_1,'.mcool.res_',ins_bin,'.window_100k.ins'))
ins_2<-fread(paste0(sample_2,'.mcool.res_',ins_bin,'.window_100k.ins'))

ins_1<-ins_1[-which(is.na(ins_1$log2_insulation_score_100000)),]
ins_2<-ins_2[-which(is.na(ins_2$log2_insulation_score_100000)),]

ins_merge<-merge(ins_1,ins_2,by=c('chrom','start'))
cor(ins_merge$log2_insulation_score_100000.x,ins_merge$log2_insulation_score_100000.y)

p_cor<-round(cor(ins_merge$log2_insulation_score_100000.x,ins_merge$log2_insulation_score_100000.y,method='pearson'),3)
s_cor<-round(cor(ins_merge$log2_insulation_score_100000.x,ins_merge$log2_insulation_score_100000.y,method='spearman'),3)

#pdf(paste0('ins_',sample_1,'_',sample_2,'_',ins_bin,'.pdf'),width=5,height=4)
p2<-ggplot(ins_merge, aes(x=ins_merge$log2_insulation_score_100000.x,y=ins_merge$log2_insulation_score_100000.y)) + geom_hex(bins =50)+ geom_abline(intercept = 0,alpha=0.2, slope = 1,linetype = "dashed") +#geom_smooth(method=lm)+ 
  scico::scale_fill_scico(palette = "lajolla",direction=-1)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold")) + xlim(-2,2) + ylim(-2,2) +
    theme(legend.direction = "horizontal",legend.position = "top", legend.box = "top")+ labs(x = NULL, y = NULL)
# labs(title=paste0("ins    x:",sample_1,", \ny:",sample_2,", hg19, bin, \npearson: ",p_cor," spearman:",s_cor))+
#p2
ins_1<-fread(paste0(sample_1,'.mcool.res_',ins_bin,'.window_100k.ins'))
ins_2<-fread(paste0(sample_2,'.mcool.res_',ins_bin,'.window_100k.ins'))
ins_1<-ins_1[ins_1$chrom=='chr7',]
ins_1$locus_id<-rownames(ins_1)
ins_2<-ins_2[ins_2$chrom=='chr7',]
ins_2$locus_id<-rownames(ins_2)
ins_1$sample<-sample_1
ins_2$sample<-sample_2
ins_line<-rbind.data.frame(ins_1,ins_2)

ins_line<-ins_line[ins_line$start>=ins_start&ins_line$end<=ins_end,]
p4<-ggplot(ins_line,aes(x=locus_id,y=log2_insulation_score_100000,group=sample,color=sample)) + geom_line(stat="identity") +
  scale_color_manual(values = c('#E64B35B2','#4DBBD5B2')) +theme_bw()+ geom_hline(yintercept = 0,linetype="dashed")+
  theme(axis.text.x=element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
  theme(legend.position="top")  + labs(x = NULL, y = NULL)+ theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

#dev.off()
dir='D:/local_project_file/hichew_2024/hicVelo_testis/pac_rou_hichew_clust/'
#dir='D:/local_project_file/hichew_2024/hichew_bulk_vs_gold_standard/'
hic1 <- import(paste0(dir,sample_1,'.mcool'),resolution=2500000)
hic2 <- import(paste0(dir,sample_2,'.mcool'),resolution=2500000)

bin=500000
hic1_chr7 <- import(CoolFile(paste0(dir,sample_1,'.mcool')), focus = 'chr7', resolution = bin)
hic2_chr7 <- import(CoolFile(paste0(dir,sample_2,'.mcool')), focus = 'chr7', resolution = bin)

scc <- get.scc(hic1["chr7"] |> as.matrix(use.scores = 'count'), 
               hic2["chr7"] |> as.matrix(use.scores = 'count'), 
               resol = 50000, h = 2, lbr = 250000, ubr = 2500000)

#pdf(paste0(sample_1,'_VS_',sample_2,'_chr7.pdf'),width=8,height=8)
p5<-plotMatrix(hic2_chr7, compare.to = hic1_chr7,caption=FALSE, use.scores = 'balanced')
#+
#  labs(title=paste0('top right: ',sample_2,' #
#
#','bottom left: ',sample_1,'
#                  
#Stratum-adjusted correlation coefficient: ',round(scc$scc[1],2)))


bin=50000
bpparam <- SerialParam(progressbar = FALSE)
hic1 <- zoom(hic1_chr7, bin) |> 
  refocus('chr7:39000000-50000000') 

hic2 <- zoom(hic2_chr7, bin) |> 
  refocus('chr7:39000000-50000000') 

hic1 <- zoom(hic1, 50000)
hic2 <- zoom(hic2, 50000)

#pdf(paste0(sample_1,'_VS_',sample_2,'_tad.pdf'),width=8,height=8)
p6<-plotMatrix(hic2, use.scores = 'balanced', 
           #borders = topologicalFeatures(hic1, "borders"),
           caption = TRUE, compare.to = hic1)
#+
#  labs(title=paste0('top right: ',sample_2,' 
#
# ','bottom left: ',sample_1))

#dev.off()
cowplot::plot_grid(p5,p1,p6,p2,nrow=1)
wrap_plots(p5,p1,p6,p2,nrow=1)
}


result <- map(sample_list, process_elements)

pdf(paste0('testis_ab_eigen_ins_tad_cor_20241126.pdf'),width=5,height=5)
cowplot::plot_grid(plotlist = result,nrow=1)
dev.off()

pdf(paste0('testis_ab_eigen_ins_tad_cor_20241126.png'),width=5,height=5)
cowplot::plot_grid(plotlist = result,nrow=1)
dev.off()

#pdf(paste0('ins_',sample_1,'_',sample_2,'_',ins_bin,'.pdf'),width=20,height=50)
#cowplot::plot_grid(p2,p4,ncol=4)
#dev.off()



#library(gridExtra)
#pdf('eigen_compare.pdf',height = 5,width = 5)
#grid.arrange(grobs=plot_list_eigen,nrow=2,ncol=3)
#dev.off()
#cowplot::plot_grid(plotlist = plot_list_eigen)

#pdf('ins_compare.pdf',height = 5,width = 5)
#grid.arrange(grobs=plot_list_ins,nrow=2,ncol=3)
#dev.off()
