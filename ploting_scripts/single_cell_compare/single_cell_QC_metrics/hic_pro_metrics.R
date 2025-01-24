library(ggplot2)
library(gridExtra)
library(ggsci)
library(reshape2)
library(purrr)
library(data.table)
library(dplyr)

mypal <- pal_npg("nrc", alpha = 0.7)(9)
color_palette<-c(mypal[1],mypal[2],mypal[3],mypal[4],mypal[5],mypal[7],mypal[8],mypal[9])

batch_f<- list(hek_sciHiC="empty_filtered_293T_sciHiC",
               hek_dipC="empty_filtered_zhichao_293tdipC_20241025",
               hek_SnHichew_200cell_b1="empty_filtered_zhichao_200cell293Tb1_20240729",
               hek_SnHichew_200cell_b2="empty_filtered_zhichao_200cell293Tb2_20240729",
               gm12878_sciHiC="empty_filtered_gm12878_sciHiC",
               gm12878_dipC="empty_filtered_dipC_gm12878_0930",
               gm12878_SnHichew_1="snHiChew_gm12878_0918",
               gm12878_SnHichew_2="empty_filtered_snHiChew_gm12878_0930")

batch_data <- map2(batch_f, names(batch_f), function(file, cell_lib) {
  dt <- fread(paste0(file, '_hic_pro_metrics.xls'))  # Read the file into a data.table
  dt[, cell_lib := cell_lib]  # Add the label column
  return(dt)})
group_df <- as.data.frame(rbindlist(batch_data))

# total reads QC
stack_df_total_reads<-group_df[,c(19,2:6)]
stack_df_total_reads_mean<-stack_df_total_reads %>% group_by(cell_lib) %>% summarize(Total_pairs_processed=mean(Total_pairs_processed))

p1<-ggplot(stack_df_total_reads_mean,aes(x=cell_lib,y=log10(Total_pairs_processed)))+
  geom_bar(stat = 'identity',alpha=0.8)+
  theme_classic()+scale_fill_manual(values='grey')+coord_flip()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

stack_df_total_reads[,c(3:6)]<-stack_df_total_reads[,c(3:6)]/stack_df_total_reads$Total_pairs_processed
stack_df_total_reads<-stack_df_total_reads[,-2]
stack_df_total_reads_mean<-stack_df_total_reads %>% group_by(cell_lib) %>% summarize(across(Unmapped_pairs:Reported_pairs, ~ mean(.x, na.rm = TRUE)))

total_reads_melted <- reshape2::melt(stack_df_total_reads_mean, id.vars = 'cell_lib')

p2<-ggplot(total_reads_melted, aes(x = cell_lib, y = value, fill = variable)) +
  geom_bar(stat = 'identity',alpha=0.8) +
  labs(x = 'Samples', y = 'mean Percentage', title = 'Total reads QC') +
  theme_classic()+scale_fill_manual(values=color_palette)+
  coord_flip()+theme(legend.position = "top")+guides(fill = guide_legend(nrow = 2))


# reported pair QC
stack_df_validpair<-group_df[,c(19,6:11)]
stack_df_validpair_mean<-stack_df_validpair %>% group_by(cell_lib) %>% summarize(Reported_pairs=mean(Reported_pairs))

p3<-ggplot(stack_df_validpair_mean,aes(x=cell_lib,y=log10(Reported_pairs)))+
  geom_bar(stat = 'identity',alpha=0.8)+
  theme_classic()+scale_fill_manual(values='grey')+coord_flip()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

stack_df_validpair[,c(3:7)]<-stack_df_validpair[,c(3:7)]/stack_df_validpair$Reported_pairs
stack_df_validpair<-stack_df_validpair[,-2]
colnames(stack_df_validpair)
stack_df_total_reads_mean<-stack_df_validpair %>% group_by(cell_lib) %>% summarize(across(Dangling_end_pairs:valid_pairs, ~ mean(.x, na.rm = TRUE)))
validpair_melted <- reshape2::melt(stack_df_total_reads_mean, id.vars = 'cell_lib')

p4<-ggplot(validpair_melted, aes(x = cell_lib, y = value, fill = variable)) +
  geom_bar(stat = 'identity',alpha=0.8) +
  labs(x = 'Samples', y = 'mean Percentage', title = 'Reported pair QC') +
  theme_classic()+scale_fill_manual(values=color_palette)+
  coord_flip()+theme(legend.position = "top")+guides(fill = guide_legend(nrow = 3))


plotlist<-list(p2,p1,p4,p3)
library(gridExtra)
pdf('single_cell_QC_metrics.pdf',width=8,height=9)
grid.arrange(grobs=plotlist,nrow=2,ncol=2,widths=c(2,1))
dev.off()

########################################################################
