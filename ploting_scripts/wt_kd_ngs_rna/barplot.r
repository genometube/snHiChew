library(ggplot2) 
library(plyr)
library(reshape2)
library(zoo)
library(pheatmap)
library(Matrix)
library(data.table)
require(magrittr)
require(readr)
require(tidyr)
require(dplyr)

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- plyr::rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df <- ToothGrowth
df$dose <- as.factor(df$dose)
head(df)
df2 <- data_summary(ToothGrowth, varname="len", 
                    groupnames=c("supp", "dose"))
# Convert dose to a factor variable
df2$dose=as.factor(df2$dose)
head(df2)
p<- ggplot(df2, aes(x=dose, y=len, fill=supp)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=len-sd, ymax=len+sd), width=.2,
                position=position_dodge(.9)) 
print(p)
p+labs(title="Tooth length per dose", x="Dose (mg)", y = "Length")+
  theme_classic() +
  scale_fill_manual(values=c('#999999','#E69F00'))

##########wt vs ko
df<-cbind.data.frame(sample=c('wt','wt','wt','siRNA','siRNA','siRNA','ko','ko','ko'),ctcf_rna_seq_exp=c(35.42,34.09,34.85,21.43,22.1,20.87,18.4,19.13,16.9))
df$sample <- as.factor(df$sample)
head(df)
df2 <- data_summary(df, varname="ctcf_rna_seq_exp", 
                    groupnames=c("sample"))
# Convert dose to a factor variable
df2$sample=as.factor(df2$sample)
head(df2)
df2$sample<-factor(df2$sample,levels=c('wt','ko','siRNA'))


p<- ggplot(df2, aes(x=sample, y=ctcf_rna_seq_exp, fill=sample)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=ctcf_rna_seq_exp, ymax=ctcf_rna_seq_exp+sd), width=.2,
                position=position_dodge(.9)) 
print(p)
pdf('RNA-seq_wt_ko_siRNA.pdf',width=4,height=7)
p+labs(title="wt vs siRNA p-value = 1.505e-05,
       wt vs ko p-value = 0.000127", x="sample", y = "RNA-seq normalized expression")+
  theme_classic() + scale_fill_manual(values=c('#999999','#E69F00','#0ABAB5'))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))
dev.off()
t.test(df[df$sample=="wt",2],df[df$sample=="siRNA",2])
t.test(df[df$sample=="wt",2],df[df$sample=="ko",2])

########## wt vs siRNA rna-seq
df<-cbind.data.frame(sample=c('wt','wt','wt','siRNA','siRNA'),ctcf_rna_seq_exp=c(35.42,34.09,34.85,11.71,11.23))
df$sample <- as.factor(df$sample)
head(df)
df2 <- data_summary(df, varname="ctcf_rna_seq_exp", 
                    groupnames=c("sample"))
# Convert dose to a factor variable
df2$sample=as.factor(df2$sample)
df2$sample<-factor(df2$sample,levels=c('wt','siRNA'))

head(df2)
p<- ggplot(df2, aes(x=sample, y=ctcf_rna_seq_exp, fill=sample)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=ctcf_rna_seq_exp, ymax=ctcf_rna_seq_exp+sd), width=.2,
                position=position_dodge(.9)) 
print(p)

pdf('RNA-seq_wt_siRNA.pdf',width=4,height=7)
p+labs(title="p-value = 1.826e-05", x="sample", y = "RNA-seq normalized expression")+
  theme_classic() + scale_fill_manual(values=c('#999999','#E69F00'))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))
t.test(df[df$sample=="wt",2],df[df$sample=="siRNA",2])
dev.off()




########## wt vs siRNA qpcr
df<-cbind.data.frame(sample=c('wt','wt','wt','siRNA','siRNA','siRNA'),ctcf_rna_seq_exp=c(1.097,0.77,1.184,0.537,0.449,0.53))
df$sample <- as.factor(df$sample)
head(df)
df2 <- data_summary(df, varname="ctcf_rna_seq_exp", 
                    groupnames=c("sample"))
# Convert dose to a factor variable
df2$sample=as.factor(df2$sample)
df2$sample<-factor(df2$sample,levels=c('wt','siRNA'))

head(df2)
p<- ggplot(df2, aes(x=sample, y=ctcf_rna_seq_exp, fill=sample)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=ctcf_rna_seq_exp, ymax=ctcf_rna_seq_exp+sd), width=.2,
                position=position_dodge(.9)) 
print(p)

pdf('qPCR_wt_siRNA.pdf',width=4,height=7)
p+labs(title="p-value = 0.02483", x="sample", y = "qPCR normalized expression")+
  theme_classic() + scale_fill_manual(values=c('#999999','#0ABAB5'))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))
t.test(df[df$sample=="wt",2],df[df$sample=="siRNA",2],alternative =c('greater'))
dev.off()

########## wt vs ko qpcr
df<-cbind.data.frame(sample=c('wt','wt','wt','ko','ko','ko'),ctcf_rna_seq_exp=c(1.14,0.991,0.887,0.74,0.839,0.746))
df$sample <- as.factor(df$sample)
head(df)
df2 <- data_summary(df, varname="ctcf_rna_seq_exp", 
                    groupnames=c("sample"))
# Convert dose to a factor variable
df2$sample=as.factor(df2$sample)
df2$sample<-factor(df2$sample,levels=c('wt','ko'))

head(df2)
p<- ggplot(df2, aes(x=sample, y=ctcf_rna_seq_exp, fill=sample)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=ctcf_rna_seq_exp, ymax=ctcf_rna_seq_exp+sd), width=.2,
                position=position_dodge(.9)) 
print(p)

pdf('qPCR_wt_ko.pdf',width=4,height=7)
p+labs(title="p-value = 0.03537", x="sample", y = "qPCR normalized expression")+
  theme_classic() + scale_fill_manual(values=c('#999999','#E69F00'))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))
t.test(df[df$sample=="wt",2],df[df$sample=="ko",2],alternative =c('greater'))
?t.test
dev.off()

#chr = "chr7"

########### line dot plot 20220413
#df<-read.table('../../../raw_table/ont_atacPoreC/hek_wt_nome3h_45M_chr7/hek_nome3h_45M_chr7_single_mol_methy_seg.xls.gz', sep='\t',header = TRUE)
#m_seg_nome3h_GpC<-df[(df$GpC_unmethy+df$GpC_methy)>4,]
#m_seg_nome3h_GpC$GpC_all<-m_seg_nome3h_GpC$GpC_methy+m_seg_nome3h_GpC$GpC_unmethy
#m_seg_nome3h_GpC
#setDT(m_seg_nome3h_GpC)
#m_seg_nome3h_GpC[,'p_val_m' := mapply(function(x,y) binom.test(x=x,n=y,p=0.02995435,alternative='greater')$p.value,GpC_methy,GpC_all)]
#m_seg_nome3h_GpC$binom_padj_m<-p.adjust(m_seg_nome3h_GpC$p_val_m,method="BH")
m_seg_wt<-fread('sca_hek_wt_chr7_single_mol_methy_seg.binom.xls.gz')
head(m_seg_wt)
read_id_all<-cbind.data.frame(read_id=m_seg_nome3h_GpC$read_id)
read_id_order<-ddply(read_id_all,.(read_id),nrow)
m_seg_nome3h_GpC<-merge(m_seg_nome3h_GpC,read_id_order,by='read_id')
head(m_seg_nome3h_GpC)
colnames(m_seg_nome3h_GpC)
colnames(m_seg_nome3h_GpC)<-c('read_id','chr','segment_start','segment_end','strand',
                              'segment_num','GpC_methy','GpC_unmethy','CpG_methy','CpG_unmethy',
                              'GpC_all','p_val_m','binom_padj_m','order')
write.table(m_seg_nome3h_GpC,file = 'm_seg_nome3h_GpC_binom.xls',quote = FALSE,sep="\t",row.names = FALSE)

m_seg_nome3h_GpC<-fread('m_seg_nome3h_GpC_binom.xls', sep='\t',header = TRUE)
# all_open / all_fragments
m_seg_nome3h_GpC<-m_seg_wt
for (i in c(1:30)) {
  order_num<-i
  all_open<-length(rownames(m_seg_nome3h_GpC[m_seg_nome3h_GpC$order==order_num&m_seg_nome3h_GpC$binom_padj_m<0.05,1]))
  all_frag<-length(rownames(m_seg_nome3h_GpC[m_seg_nome3h_GpC$order==order_num,1]))
  print(all_open/all_frag)
}

head(m_seg_nome3h_GpC)
id_open<-m_seg_nome3h_GpC[m_seg_nome3h_GpC$binom_padj_m<0.05,1]
id_open<-as.data.frame(id_open)
colnames(id_open)<-'read_id'
id_open_freq<-ddply(id_open,.(read_id),nrow)
id_open_freq_2<-cbind.data.frame(read_id=id_open_freq[id_open_freq$V1>=2,1])
m_seg_nome3h_GpC_open<-merge(id_open_freq_2,m_seg_nome3h_GpC,by='read_id')
order_to_2open_ratio<-cbind.data.frame(order=1,ratio=0,sample_size=0)
for (i in c(2:18)) {
  order_num<-i
  all_open<-length(m_seg_nome3h_GpC_open[m_seg_nome3h_GpC_open$order==order_num&m_seg_nome3h_GpC_open$binom_padj_m<0.05,1])
  all_frag<-length(rownames(m_seg_nome3h_GpC[m_seg_nome3h_GpC$order==order_num,1]))
  print(all_open/all_frag)
  tmp <-cbind.data.frame(order=i,ratio=all_open/all_frag,sample_size=all_frag)
  order_to_2open_ratio<-rbind.data.frame(order_to_2open_ratio,tmp)
}
order_to_2open_ratio<-order_to_2open_ratio[-1,]
cor(order_to_2open_ratio[,1],order_to_2open_ratio[,2])
loess(order_to_2open_ratio$order ~ order_to_2open_ratio$ratio,span=1,degree=2)

order_to_2open_ratio.lo <- loess(order_to_2open_ratio$order ~ order_to_2open_ratio$ratio,span=1,degree=2)
summary(loess(order ~ ratio, data=order_to_2open_ratio, span=1,model=TRUE) )

library(mgcv)
summary(gam(formula = ratio ~ s(order,k=3), data=order_to_2open_ratio))

pdf('wt_order_to_2open_ratio.pdf') 
ggplot(order_to_2open_ratio,aes(x=order,y=ratio)) + geom_point(color='#FF4500') + ylim(0,0.15)+theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
  stat_smooth(method = "gam",se = FALSE, formula = y ~ s(x,k=3), size=1,color=c('#FF4500')) +
  scale_linetype_discrete(name="type", breaks=c("eigenvector", "insulation_score"))+
  labs(title="SCA_wt, GAM_model_R_square=0.984
       x-axis: order,
       y-axis: open concatemer ratio (>= 2 open chromatin per concatemer)") +
  scale_x_continuous(breaks=c(2:18))
dev.off()
#geom_smooth(method = "loess",se = FALSE,span = 1, size=0.5,color=c('#FF4500')) +
# ,colour = type,linetype=type +scale_color_manual(values=c("#A52A2A","#FF4500"))

sum_sample_size<-sum(order_to_2open_ratio$sample_size)
pdf('wt_order_to_concatemer_pct.pdf') 
ggplot(order_to_2open_ratio,aes(x=order,y=sample_size/sum_sample_size)) + geom_point(color='#FF4500') + ylim(0,0.3)+theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
  geom_smooth(method = "loess",se = FALSE,span = 0.5, size=1,color=c('#FF4500')) +
  scale_linetype_discrete(name="type", breaks=c("eigenvector", "insulation_score"))+
  labs(title="SCA_wt,
       x-axis: order,
       y-axis: concatemer pct") +
scale_x_continuous(breaks=c(2:18))
dev.off()

order_to_2open_ratio_ko<-fread('ko_order_to_2open_ratio.xls', sep='\t',header = TRUE)
order_to_2open_ratio_ko<-order_to_2open_ratio_ko[-1,]
order_to_2open_ratio$group<-'SCA_wt'
order_to_2open_ratio$concatemer_pct<-order_to_2open_ratio$sample_size/sum(order_to_2open_ratio$sample_size)
order_to_2open_ratio_ko$group<-'SCA_ctcf_ko'
order_to_2open_ratio_ko$concatemer_pct<-order_to_2open_ratio_ko$sample_size/sum(order_to_2open_ratio_ko$sample_size)
order_to_2open_ratio_ko_wt_ko<-rbind.data.frame(order_to_2open_ratio,order_to_2open_ratio_ko)

pdf('wt_vs_ko_order_to_2open_ratio.pdf') 
ggplot(order_to_2open_ratio_ko_wt_ko,aes(x=order,y=ratio,color = group,linetype=group)) + geom_point(aes(color=group)) + ylim(0,0.25)+theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
  stat_smooth(method = "gam",se = FALSE, formula = y ~ s(x,k=3), size=1) +
  scale_color_manual(values=c("#DAA520","#FF4500"))+scale_linetype_manual(values=c("dashed", "solid"))+

  labs(title="SCA_wt_vs_ctcfko, GAM_model fit
       x-axis: order,
       y-axis: open concatemer ratio (>= 2 open chromatin per concatemer)") +
  scale_x_continuous(breaks=c(2:18))
dev.off()
#  geom_smooth(method = "loess",se = FALSE,span = 1, size=0.7) +   scale_linetype_discrete(name="type")+

pdf('wt_vs_ko_order_to_concatemer_pct.pdf') 
ggplot(order_to_2open_ratio_ko_wt_ko,aes(x=order,y=concatemer_pct,color=group,linetype=group)) + geom_point(aes(color=group)) + ylim(0,0.3)+theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
  geom_smooth(method = "loess",se = FALSE,span = 0.5, size=1) +
  scale_color_manual(values=c("#DAA520","#FF4500"))+scale_linetype_manual(values=c("dashed", "solid"))+
  labs(title="SCA_wt,
       x-axis: order,
       y-axis: concatemer pct") +
  scale_x_continuous(breaks=c(2:18))
dev.off()

for (i in c(1:5)) {
  order_num_min<-2^i
  order_num_max<-2^(i+1)
  all_open<-length(m_seg_nome3h_GpC_open[m_seg_nome3h_GpC_open$order>=order_num_min&m_seg_nome3h_GpC_open$order<=order_num_max&m_seg_nome3h_GpC_open$binom_padj_m<0.05,1])
  all_frag<-length(m_seg_nome3h_GpC[m_seg_nome3h_GpC$order>=order_num_min&m_seg_nome3h_GpC$order<=order_num_max,1])
  print(all_open/all_frag)
}

order_num<-14
all_open<-length(m_seg_nome3h_GpC[m_seg_nome3h_GpC$order==order_num&m_seg_nome3h_GpC$binom_padj_m<0.05,1])
all_frag<-length(m_seg_nome3h_GpC[m_seg_nome3h_GpC$order==order_num,1])
all_open/all_frag

order_num_min<-2
order_num_max<-30
all_open<-length(m_seg_nome3h_GpC[m_seg_nome3h_GpC$order>=order_num_min&m_seg_nome3h_GpC$order<=order_num_max&m_seg_nome3h_GpC$binom_padj_m<0.05,1])
all_frag<-length(m_seg_nome3h_GpC[m_seg_nome3h_GpC$order>=order_num_min&m_seg_nome3h_GpC$order<=order_num_max,1])
all_open/all_frag

m_seg_nome3h_GpC<-fread('m_seg_nome3h_GpC_binom.xls', sep='\t',header = TRUE)
length(which(m_seg_nome3h_GpC$binom_padj_m<0.05))
#df_chr = df[df$V1==chr,c(2,3)]
seg_pos_GpC_ratio<-cbind.data.frame(read_id=m_seg_nome3h_GpC$read_id,
                                    pos=(m_seg_nome3h_GpC$segment_start+m_seg_nome3h_GpC$segment_end)/2,
                                    GpC_ratio=m_seg_nome3h_GpC$GpC_methy/(0.000001+m_seg_nome3h_GpC$GpC_methy+m_seg_nome3h_GpC$GpC_unmethy),
                                    padj=m_seg_nome3h_GpC$binom_padj_m)
seg_bin_GpC_ratio<-seg_pos_GpC_ratio


########### line dot plot 20220414
# selected region
bin=1000
locus_1<-75900000
locus_2<-76000000
seg_bin_GpC_ratio<-seg_pos_GpC_ratio[seg_pos_GpC_ratio$pos>locus_1&seg_pos_GpC_ratio$pos<locus_2,]
seg_bin_GpC_ratio$pos <- floor(seg_bin_GpC_ratio$pos/bin)
min(seg_bin_GpC_ratio[seg_bin_GpC_ratio$padj<0.05,3])
seg_bin_GpC_ratio[seg_bin_GpC_ratio$padj>=0.05,3]<- -0.3
  
m_seg_read_id_order<-ddply(seg_bin_GpC_ratio,.(read_id),nrow)
head(m_seg_read_id_order)
#keep order > 1 and order < 5 read_id
#m_seg_read_id_order_2_<-m_seg_read_id_order[m_seg_read_id_order$V1>1&m_seg_read_id_order$V1<20,]
m_seg_read_id_order_2_<-m_seg_read_id_order[m_seg_read_id_order$V1>1,]
seg_bin_GpC_ratio_order_2_<-merge(m_seg_read_id_order_2_,seg_bin_GpC_ratio,by='read_id')
seg_bin_GpC_ratio_order_2_sort <- seg_bin_GpC_ratio_order_2_
seg_bin_GpC_ratio_order_2_sort_pos_mean<-aggregate(pos~read_id,data=seg_bin_GpC_ratio_order_2_sort , mean)
seg_bin_GpC_ratio_order_2_sort<-merge(seg_bin_GpC_ratio_order_2_sort,seg_bin_GpC_ratio_order_2_sort_pos_mean,by='read_id')
colnames(seg_bin_GpC_ratio_order_2_sort)<-c('read_id','order','pos','GpC_ratio','padj','pos_mean')
# sort by order and pos mean
seg_bin_GpC_ratio_order_2_sort <- seg_bin_GpC_ratio_order_2_sort[order(seg_bin_GpC_ratio_order_2_sort$order,seg_bin_GpC_ratio_order_2_sort$pos_mean),]

seg_bin_GpC_ratio_order_2_sort$read_id<-factor(seg_bin_GpC_ratio_order_2_sort$read_id
                                               ,levels=seg_bin_GpC_ratio_order_2_sort[!duplicated(seg_bin_GpC_ratio_order_2_sort$read_id),1])
seg_bin_GpC_ratio_order_2_sort$order<-as.factor(seg_bin_GpC_ratio_order_2_sort$order)
head(seg_bin_GpC_ratio_order_2_sort)

#seg_bin_GpC_ratio_order_2_sort$order<-as.character(seg_bin_GpC_ratio_order_2_sort$order)
seg_bin_GpC_ratio_order_2_sort[seg_bin_GpC_ratio_order_2_sort$GpC_ratio>0.4,4] <- 0.3

#pdf(paste0("SCA_WT_chr7_",locus_1,'-',locus_2,'_bin_',bin,'_f_',feature_start,'-',feature_end,'.pdf')) # size=GpC_ratio,
ggplot(seg_bin_GpC_ratio_order_2_sort,aes(x=pos,y=read_id)) + geom_point(alpha=0.9,aes(shape = order,color = GpC_ratio,size=GpC_ratio)) +
  scale_shape_manual(values=1:nlevels(seg_bin_GpC_ratio_order_2_sort$order))+
  scale_size_continuous(range = c(1,3))+scale_color_continuous(c(-0.2,0.2))+theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "grey"),axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  scale_color_gradient2(low = "grey30", mid = "white", high = "red")+
  geom_line(stat="smooth",method = "lm",color='darkgrey',alpha = 0.05)+
  labs(title=paste0("SCA_WT, chr7:",locus_1,'-',locus_2,'   bin_size: ',bin))
  #geom_vline(xintercept  = c(26237883/1000-0.5,26242631/1000), col = "darkgreen",linetype = "dashed")+
  #geom_vline(xintercept  = c(26190000/1000,26192000/1000,), col = "darkblue",linetype = "dashed")+
dev.off() 

########### line dot plot 20220414
############################ selected region with 1 feature
locus_1<-82900000 
locus_2<-85000000 
feature_start<-83277660 
feature_end<-83278740
bin=1000
locus_1<-139000000
locus_2<-141100000
feature_start<-140396000
feature_end<-140398000
feature_id<-cbind.data.frame(read_id=seg_pos_GpC_ratio[seg_pos_GpC_ratio$pos>feature_start&seg_pos_GpC_ratio$pos<feature_end,1])

seg_bin_GpC_ratio<-seg_pos_GpC_ratio[seg_pos_GpC_ratio$pos>locus_1&seg_pos_GpC_ratio$pos<locus_2,]
seg_bin_GpC_ratio$pos <- floor(seg_bin_GpC_ratio$pos/bin)
min(seg_bin_GpC_ratio[seg_bin_GpC_ratio$padj<0.05,3])
seg_bin_GpC_ratio[seg_bin_GpC_ratio$padj>=0.05,3]<- -0.3
seg_bin_GpC_ratio<-merge(feature_id,seg_bin_GpC_ratio,by='read_id')
length(feature_id$read_id)
head(feature_id)
head(seg_bin_GpC_ratio)

m_seg_read_id_order<-ddply(seg_bin_GpC_ratio,.(read_id),nrow)
head(m_seg_read_id_order)
#keep order > 1 and order < 5 read_id
#m_seg_read_id_order_2_<-m_seg_read_id_order[m_seg_read_id_order$V1>1&m_seg_read_id_order$V1<5,]
m_seg_read_id_order_2_<-m_seg_read_id_order[m_seg_read_id_order$V1>1,]
seg_bin_GpC_ratio_order_2_<-merge(m_seg_read_id_order_2_,seg_bin_GpC_ratio,by='read_id')
seg_bin_GpC_ratio_order_2_sort <- seg_bin_GpC_ratio_order_2_
seg_bin_GpC_ratio_order_2_sort_pos_mean<-aggregate(pos~read_id,data=seg_bin_GpC_ratio_order_2_sort , mean)
seg_bin_GpC_ratio_order_2_sort<-merge(seg_bin_GpC_ratio_order_2_sort,seg_bin_GpC_ratio_order_2_sort_pos_mean,by='read_id')
colnames(seg_bin_GpC_ratio_order_2_sort)<-c('read_id','order','pos','GpC_ratio','padj','pos_mean')
#remove open <= 1 read_id
head(seg_bin_GpC_ratio_order_2_sort)
id_open_freq<-as.data.frame(table(seg_bin_GpC_ratio_order_2_sort[seg_bin_GpC_ratio_order_2_sort$padj<0.05,1]))
#id_open_freq_2_<-id_open_freq[id_open_freq$Freq>1,]
id_open_freq_2_<-id_open_freq
id_open_2_<-cbind.data.frame(read_id=id_open_freq_2_$Var1)
#filter_id<-seg_bin_GpC_ratio_order_2_sort[seg_bin_GpC_ratio_order_2_sort$padj<0.05,1]
#filter_id<-filter_id[!duplicated(filter_id)]
#filter_id<-cbind.data.frame(read_id=filter_id)
seg_bin_GpC_ratio_order_2_sort<-merge(id_open_2_,seg_bin_GpC_ratio_order_2_sort,by='read_id')
#remove feature closed read_id
feature_open_id<-cbind.data.frame(read_id=seg_bin_GpC_ratio_order_2_sort[seg_bin_GpC_ratio_order_2_sort$pos>=feature_start/bin&seg_bin_GpC_ratio_order_2_sort$pos<=feature_end/bin&seg_bin_GpC_ratio_order_2_sort$padj<0.05,1])
seg_bin_GpC_ratio_order_2_sort<-merge(feature_open_id,seg_bin_GpC_ratio_order_2_sort,by='read_id')
# sort by order and pos mean
seg_bin_GpC_ratio_order_2_sort <- seg_bin_GpC_ratio_order_2_sort[order(seg_bin_GpC_ratio_order_2_sort$order,seg_bin_GpC_ratio_order_2_sort$pos_mean),]

seg_bin_GpC_ratio_order_2_sort$read_id<-factor(seg_bin_GpC_ratio_order_2_sort$read_id
  ,levels=seg_bin_GpC_ratio_order_2_sort[!duplicated(seg_bin_GpC_ratio_order_2_sort$read_id),1])
seg_bin_GpC_ratio_order_2_sort$order<-as.factor(seg_bin_GpC_ratio_order_2_sort$order)
head(seg_bin_GpC_ratio_order_2_sort)
open_frag_num<-length(rownames(seg_bin_GpC_ratio_order_2_sort[seg_bin_GpC_ratio_order_2_sort$padj<0.5,]))
total_frag_num<-length(rownames(seg_bin_GpC_ratio_order_2_sort))

#seg_bin_GpC_ratio_order_2_sort$order<-as.character(seg_bin_GpC_ratio_order_2_sort$order)
seg_bin_GpC_ratio_order_2_sort[seg_bin_GpC_ratio_order_2_sort$GpC_ratio>0,4] <- 0.25
seg_bin_GpC_ratio_order_2_sort$GpC_ratio<-as.factor(seg_bin_GpC_ratio_order_2_sort$GpC_ratio)
#scale_color_continuous(c(-0.2,0.2))+scale_color_gradient2(low = "grey30", mid = "white", high = "red")+scale_size_continuous(range = c(3,4))+
pdf(paste0("whole_region_SCA_WT_chr7_",locus_1,'-',locus_2,'_bin_',bin,'_f_',feature_start,'-',feature_end,'.pdf')) # size=GpC_ratio,
ggplot(seg_bin_GpC_ratio_order_2_sort,aes(x=pos,y=read_id)) + geom_point(alpha=0.8,aes(shape = order,color = GpC_ratio),stroke=1) +
  scale_shape_manual(values=1:nlevels(seg_bin_GpC_ratio_order_2_sort$order))+
  theme_bw()+scale_color_manual(values=c('grey30','red'))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "lightgrey"),axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  geom_line(stat="smooth",method = "lm",color='grey',alpha = 0.8)+xlim(locus_1/1000,locus_2/1000)+
  labs(title=paste0("SCA_WT chr7:",locus_1,'-',locus_2,' bin_size: ',bin,' feature:',feature_start,'-',feature_end),
       subtitle = paste0('open_num: ',open_frag_num,'    total_num: ',total_frag_num))
  #geom_vline(xintercept  = c(feature_start/1000-0.5,feature_end/1000), col = "darkgreen",linetype = "dashed")
dev.off() 

open_id<-as.data.frame(seg_bin_GpC_ratio_order_2_sort[seg_bin_GpC_ratio_order_2_sort$pos>=feature_start/bin&seg_bin_GpC_ratio_order_2_sort$pos<=feature_end/bin&seg_bin_GpC_ratio_order_2_sort$padj<0.05,1])
colnames(open_id)<-'read_id'
feature_open<-merge(open_id,seg_bin_GpC_ratio_order_2_sort,by='read_id')
freq<-as.data.frame(table(feature_open[feature_open$padj<0.05,1]))
length(which(freq$Freq >1))

########### line dot plot 20220414
# selected enhancer promoter
bin=1000
locus_1<-139700435
locus_2<-141205680
feature_start<-104614386
feature_end<-104656455
feature_start_2<-103840000
feature_end_2<-103880000
# enhancer 42312800 42325380
# enhancer 42312820 42325201
# enhancer 42784940 42785320
# promoter 42951060 42952180
# promoter 42971560 42972180
#feature_id<-cbind.data.frame(read_id=seg_pos_GpC_ratio[seg_pos_GpC_ratio$pos>feature_start&seg_pos_GpC_ratio$pos<feature_end,1])
feature_id_1<-cbind.data.frame(read_id=seg_pos_GpC_ratio[seg_pos_GpC_ratio$pos>feature_start&seg_pos_GpC_ratio$pos<feature_end,1])
feature_id_2<-cbind.data.frame(read_id=seg_pos_GpC_ratio[seg_pos_GpC_ratio$pos>feature_start_2&seg_pos_GpC_ratio$pos<feature_end_2,1])
feature_id_merge<-merge(feature_id_1,feature_id_2,by='read_id')
seg_bin_GpC_ratio<-seg_pos_GpC_ratio[seg_pos_GpC_ratio$pos>locus_1&seg_pos_GpC_ratio$pos<locus_2,]
seg_bin_GpC_ratio$pos <- floor(seg_bin_GpC_ratio$pos/bin)
min(seg_bin_GpC_ratio[seg_bin_GpC_ratio$padj<0.05,3])
seg_bin_GpC_ratio[seg_bin_GpC_ratio$padj>=0.05,3]<- -0.3
seg_bin_GpC_ratio<-merge(feature_id_merge,seg_bin_GpC_ratio,by='read_id')

head(feature_id)
head(seg_bin_GpC_ratio)
m_seg_read_id_order<-ddply(seg_bin_GpC_ratio,.(read_id),nrow)
head(m_seg_read_id_order)
#keep order > 1 and order < 5 read_id
#m_seg_read_id_order_2_<-m_seg_read_id_order[m_seg_read_id_order$V1>1&m_seg_read_id_order$V1<6,]
m_seg_read_id_order_2_<-m_seg_read_id_order[m_seg_read_id_order$V1>1,]

seg_bin_GpC_ratio_order_2_sort <- seg_bin_GpC_ratio_order_2_
seg_bin_GpC_ratio_order_2_<-merge(m_seg_read_id_order_2_,seg_bin_GpC_ratio,by='read_id')
seg_bin_GpC_ratio_order_2_sort_pos_mean<-aggregate(pos~read_id,data=seg_bin_GpC_ratio_order_2_sort , mean)
seg_bin_GpC_ratio_order_2_sort<-merge(seg_bin_GpC_ratio_order_2_sort,seg_bin_GpC_ratio_order_2_sort_pos_mean,by='read_id')
colnames(seg_bin_GpC_ratio_order_2_sort)<-c('read_id','order','pos','GpC_ratio','padj','pos_mean')
# sort by order and pos mean
seg_bin_GpC_ratio_order_2_sort <- seg_bin_GpC_ratio_order_2_sort[order(seg_bin_GpC_ratio_order_2_sort$order,seg_bin_GpC_ratio_order_2_sort$pos_mean),]

seg_bin_GpC_ratio_order_2_sort$read_id<-factor(seg_bin_GpC_ratio_order_2_sort$read_id
                                               ,levels=seg_bin_GpC_ratio_order_2_sort[!duplicated(seg_bin_GpC_ratio_order_2_sort$read_id),1])
seg_bin_GpC_ratio_order_2_sort$order<-as.factor(seg_bin_GpC_ratio_order_2_sort$order)
head(seg_bin_GpC_ratio_order_2_sort)

#seg_bin_GpC_ratio_order_2_sort$order<-as.character(seg_bin_GpC_ratio_order_2_sort$order)
seg_bin_GpC_ratio_order_2_sort[seg_bin_GpC_ratio_order_2_sort$GpC_ratio>0,4] <- 0.25
pdf(paste0("two_feature_SCA_WT_chr7_",locus_1,'-',locus_2,'_bin_',bin,'_f_',feature_start,'-',feature_end,'.pdf')) # size=GpC_ratio,
ggplot(seg_bin_GpC_ratio_order_2_sort,aes(x=pos,y=read_id)) + geom_point(alpha=0.7,aes(shape = order,color = GpC_ratio)) +
  scale_shape_manual(values=1:nlevels(seg_bin_GpC_ratio_order_2_sort$order))+
  scale_size_continuous(range = c(1,3))+scale_color_continuous(c(-0.2,0.2))+theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "lightgrey"),axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  scale_color_gradient2(low = "grey30", mid = "white", high = "red")+
  geom_line(stat="smooth",method = "lm",color='lightgrey',alpha = 0.1)+
  labs(title=paste0("two_feature_SCA_WT, chr7:",locus_1,'-',locus_2,'   bin_size: ',bin,' feature:',feature_start,'-',feature_end))
dev.off() 


########### line dot plot 20220414



dat_mtx<-dcast(seg_bin_GpC_ratio_order_2_,read_id~pos,value.var = "GpC_ratio",fun.aggregate = sum)
dat_mtx_df <- as.data.frame(dat_mtx)
dat_mtx_df2 <- dat_mtx_df[,-1]
row.names(dat_mtx_df2)<- dat_mtx_df$V2
#dat_mtx_df2[dat_mtx_df2>0.03]<-0.03
#dat_mtx_df2_log2 <-log10(dat_mtx_df2+1)*10
#dat_mtx_df2_log2[dat_mtx_df2_log2>0.3]<-0.3
pheatmap(color=colorRampPalette(c("white","red"))(50),dat_mtx_df2,border_color=NA,cluster_rows=T,cluster_cols=F,
         show_rownames=F,show_colnames=F,main = paste0('chr7:',locus_1,'-',locus_2,'   bin_size: ',bin))
# ,file=out_file
#pheatmap(color=colorRampPalette(c("#F7F7F7", 'yellow','red',"darkred"))(50),dat_mtx_df2,border_color=NA,cluster_rows=F,cluster_cols=F,
#         show_rownames=F,show_colnames=F,main = paste0(chr,':',locus_1,'-',locus_2,'   bin_size: ',bin,'    methy/all: bin1*bin2/(bin1+bin2)'))




