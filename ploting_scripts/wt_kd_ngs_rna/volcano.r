library(ggplot2)
library(plyr)
library(reshape2)
library(zoo)
library(Matrix)
library(data.table)
library(dplyr)

wt_kd_RNA<-read.table('groupA-vs-groupC.gene.DEseq2.Filter.xls',header = T)
wt_kd_RNA_subset<-wt_kd_RNA[wt_kd_RNA$GeneID=='10664',c(20,21,22,36,37,38)]
df<-reshape2::melt(wt_kd_RNA_subset)

df$condition<-sub(".*_(.).*", "\\1", df$variable)
colnames(df)<-c('sample','expression','condition')

#df$sample <- as.factor(df$sample)

summary_data <- df %>%
  group_by(condition) %>%
  summarise(
    mean_expression = mean(expression),
    se_expression = sd(expression) / sqrt(n())
  )

t_test_result <- t.test(expression ~ condition, data = df)

pdf('wt_kd_RNA_ctcf_bar.pdf',width = 6,height=8)
# Create bar plot with error bars
ggplot(summary_data, aes(x = condition, y = mean_expression, fill = condition)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.5,alpha=0.5) +
  geom_errorbar(aes(ymin = mean_expression - se_expression, ymax = mean_expression + se_expression),
                width = 0.2, position = position_dodge(0.7)) +
  labs(title = "N=3, A:hek WT, C:hek CTCF KD, Gene Expression Comparison", y = "Mean Expression", x = "Condition") +
  theme_classic()  + scale_fill_manual(values=c('#999999','#0ABAB5'))+
  annotate("text", x = 2, y = max(summary_data$mean_expression) + 0.1,
           label = paste("p-value:", format(t_test_result$p.value, digits = 3)))
dev.off()

summary_df <- df %>%
  group_by(condition) %>%
  summarise(mean_value = mean(expression),
            se = sd(expression) / sqrt(n()))
mypal <- pal_npg("nrc", alpha = 0.7)(9)

p2 <- ggplot(df, aes(x = condition, y = expression,color=condition)) +
  geom_point(position = position_jitter(width = 0.1), size = 3) +  # Add dots with jitter
  geom_errorbar(data = summary_df, aes(y = mean_value, ymin = mean_value - se, ymax = mean_value + se), width = 0.2) +  # Add error bars
  labs(title="valid / report ratio")+
  theme_classic()+theme(axis.text.x = element_text(angle = 0, hjust = 1))+ scale_color_manual(values=c('#999999','#0ABAB5'))
p2

##########wt vs ko
#df<-cbind.data.frame(sample=c('wt','wt','wt','siRNA','siRNA','siRNA','ko','ko','ko'),ctcf_rna_seq_exp=c(35.42,34.09,34.85,21.43,22.1,20.87,18.4,19.13,16.9))
df$sample <- as.factor(df$sample)
typeof(as.numeric(df$sample))
df$sample
head(df)
df2 <- data_summary(df, varname="ctcf_rna_seq_exp", 
                    groupnames=c("sample"))
# Convert dose to a factor variable
df2$sample=as.factor(df2$sample)
head(df2)
df2$sample<-factor(df2$sample,levels=c('wt','kd'))


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

wt_kd_RNA<-read.table('groupA-vs-groupC.gene.DEseq2.Filter.xls',header = T)
wt_kd_RNA_key<-wt_kd_RNA[,c(1,9:12)]
################################## 20240902
ggplot(wt_kd_RNA,aes(x=wt_kd_RNA$Software_Fold_Change,y=-log10(wt_kd_RNA$Corrected_Significant_Difference_Value)))+
  geom_point(alpha=0.1)+xlim(-5,5)+theme_classic()+
  geom_hline(yintercept=c(1.30103),linetype="dotted")+geom_vline(xintercept=c(1,-1),linetype="dotted")
wt_ko_sig<-rbind.data.frame(cbind.data.frame(group='sig_ko_wt_up',
                                             log10_padj=wt_kd_RNA[wt_kd_RNA$Software_Fold_Change>=1&-log10(wt_kd_RNA$Corrected_Significant_Difference_Value)>1.3,11],
                                             fc_ko_over_wt =wt_kd_RNA[wt_kd_RNA$Software_Fold_Change>=1&-log10(wt_kd_RNA$Corrected_Significant_Difference_Value)>1.3,12]),
                            cbind.data.frame(group='sig_ko_wt_down',
                                             log10_padj=wt_kd_RNA[wt_kd_RNA$Software_Fold_Change<= -1&-log10(wt_kd_RNA$Corrected_Significant_Difference_Value)>1.3,11],
                                             fc_ko_over_wt =wt_kd_RNA[wt_kd_RNA$Software_Fold_Change<= -1&-log10(wt_kd_RNA$Corrected_Significant_Difference_Value)>1.3,12]),
                            cbind.data.frame(group='wt_ko_no_sig',
                                             log10_padj=wt_kd_RNA[abs(wt_kd_RNA$Software_Fold_Change)<1|-log10(wt_kd_RNA$Corrected_Significant_Difference_Value)<1.3,11],
                                             fc_ko_over_wt =wt_kd_RNA[abs(wt_kd_RNA$Software_Fold_Change)<1|-log10(wt_kd_RNA$Corrected_Significant_Difference_Value)<1.3,12]))
wt_ko_sig$log10_padj<--log10(wt_ko_sig$log10_padj)
table(wt_ko_sig$group)
wt_ko_sig$group<-factor(wt_ko_sig$group,levels=c('wt_ko_no_sig','sig_ko_wt_down','sig_ko_wt_up'))

pdf('wt_ko_fisher_fc.pdf')
p1<-ggplot(wt_ko_sig, aes(x=fc_ko_over_wt,y=log10_padj,color=group)) + geom_point(alpha=0.3)+ 
  scale_color_manual(values=c("grey","#0ABAB5",'#673147'))+geom_hline(yintercept=c(1.30103),linetype="dotted")+geom_vline(xintercept=c(1,-1),linetype="dotted")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold")) + xlim(-4,4) + ylim(0,15) +
  labs(title="hg19, no_sig: 676, sig_ko_wt_up: 571, sig_ko_wt_down: 114")
dev.off()

svg('wt_ko_fisher_fc.svg')
p1
dev.off()
################################## 20240902

dat <- data.frame("smoke_no" = c(7, 0),"smoke_yes" = c(2, 5),row.names = c("Athlete", "Non-athlete"),stringsAsFactors = FALSE)
colnames(dat) <- c("Non-smoker", "Smoker")
fisher.test(data.frame("smoke_no" = c(7, 0),"smoke_yes" = c(2, 5),row.names = c("Athlete", "Non-athlete"),stringsAsFactors = FALSE))[1]
dat
test_df<-cbind.data.frame(sample='s1',wt_m=1,wt_um=10,ko_m=5,ko_um=2)
test_df
fisher.test(data.frame("wt" = c(test_df$wt_m, test_df$wt_um),"ko" = c(test_df$ko_m, test_df$ko_um),row.names = c("m", "um"),stringsAsFactors = FALSE))$p.value
data.frame("wt" = c(test_df$wt_m, test_df$wt_um),"ko" = c(test_df$ko_m, test_df$ko_um),row.names = c("m", "um"),stringsAsFactors = FALSE)
setDT(test_df)
test_df[,'p_val_fisher' := mapply(function(w,x,y,z) fisher.test(data.frame("wt" = c(w, x),"ko" = c(y, z),row.names = c("m","um"),stringsAsFactors = FALSE))$p.value,wt_m,wt_um,ko_m,ko_um)]
test_df$padj_fisher<-p.adjust(test_df$padj_fisher,method="BH")


wt_uniq<-fread('wt_unique_gpc_m_all_gpc.bedgraph')
wt_uniq_ko<-fread('wt_unique_ctcfko_hg19mt_gpc_m_all_gpc.bedgraph')
ko_uniq<-fread('ko_unique_gpc_m_all_gpc.bedgraph')
ko_uniq_wt<-fread('ko_unique_hekNome3h_45M_hg19mt_gpc_m_all_gpc.bedgraph')
wt_ko_1<-merge(wt_uniq,wt_uniq_ko,by=c('V1','V2','V3'))
wt_ko_2<-merge(ko_uniq_wt,ko_uniq,by=c('V1','V2','V3'))
wt_ko<-rbind.data.frame(wt_ko_1,wt_ko_2)
colnames(wt_ko)<-c('chr','start','end','wt_m','wt_um','ko_m','ko_um')
setDT(wt_ko)
wt_ko[,'p_val_fisher' := mapply(function(w,x,y,z) fisher.test(data.frame("wt" = c(w, x),"ko" = c(y, z),row.names = c("m","um"),stringsAsFactors = FALSE))$p.value,wt_m,wt_um,ko_m,ko_um)]
wt_ko$padj_fisher<-p.adjust(wt_ko$p_val_fisher,method="BH")
wt_ko$fc_ko_over_wt<-log2((wt_ko$ko_m/(wt_ko$ko_m+wt_ko$ko_um+0.00001))/(wt_ko$wt_m/(wt_ko$wt_m+wt_ko$wt_um+0.00001)))
wt_ko$log10_padj_fisher<- -log10(wt_ko$padj_fisher)
ggplot(wt_ko,aes(x=fc_ko_over_wt,y=log10_padj_fisher))+geom_point(alpha=0.1)+xlim(-5,5)+
  geom_hline(yintercept=c(1.30103),linetype="dotted")+geom_vline(xintercept=c(1,-1),linetype="dotted")
head(wt_ko)
write.table(wt_ko,file = 'wt_ko_fisher_exact.xls',quote = FALSE,sep="\t",row.names = FALSE)
write.table(wt_ko[abs(wt_ko$fc_ko_over_wt>=1)&wt_ko$log10_padj_fisher>1,],file = 'wt_ko_fisher_exact_sig.xls',quote = FALSE,sep="\t",row.names = FALSE)
length(rownames(wt_ko[abs(wt_ko$fc_ko_over_wt>=1)&wt_ko$log10_padj_fisher>1,]))
################################## 20220927
wt_ko_sig<-rbind.data.frame(cbind.data.frame(group='sig_ko_wt_up',
                                             log10_padj_fisher=wt_ko[wt_ko$fc_ko_over_wt>=1&wt_ko$log10_padj_fisher>1.3,11],
                                             fc_ko_over_wt =wt_ko[wt_ko$fc_ko_over_wt>=1&wt_ko$log10_padj_fisher>1.3,10]),
                            cbind.data.frame(group='sig_ko_wt_down',
                                             log10_padj_fisher=wt_ko[wt_ko$fc_ko_over_wt<= -1&wt_ko$log10_padj_fisher>1.3,11],
                                             fc_ko_over_wt =wt_ko[wt_ko$fc_ko_over_wt<= -1&wt_ko$log10_padj_fisher>1.3,10]),
                            cbind.data.frame(group='wt_ko_no_sig',
                                             log10_padj_fisher=wt_ko[abs(wt_ko$fc_ko_over_wt)<1|wt_ko$log10_padj_fisher<1.3,11],
                                             fc_ko_over_wt =wt_ko[abs(wt_ko$fc_ko_over_wt)<1|wt_ko$log10_padj_fisher<1.3,10]))
wt_ko_sig$group<-factor(wt_ko_sig$group,levels=c('wt_ko_no_sig','sig_ko_wt_down','sig_ko_wt_up'))
pdf('wt_ko_fisher_fc.pdf')
ggplot(wt_ko_sig, aes(x=fc_ko_over_wt,y=log10_padj_fisher,color=group)) + geom_point(alpha=0.15)+ 
  scale_color_manual(values=c("grey","#0ABAB5",'#673147'))+geom_hline(yintercept=c(1.30103),linetype="dotted")+geom_vline(xintercept=c(1,-1),linetype="dotted")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold")) + xlim(-4,4) + ylim(0,15) +
  labs(title="hg19, no_sig: 37256, sig_ko_wt_up: 3457, sig_ko_wt_down: 932")
dev.off()
length(which(wt_ko_sig$group=='wt_ko_no_sig'))
length(which(wt_ko_sig$group=='sig_ko_wt_down'))
length(which(wt_ko_sig$group=='sig_ko_wt_up'))

################################## 20220928
