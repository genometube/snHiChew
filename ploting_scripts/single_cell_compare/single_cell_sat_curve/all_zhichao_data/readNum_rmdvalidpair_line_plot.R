library(data.table)
library(ggplot2)
library(purrr)
library(dplyr)
library(ggsci)
mypal <- pal_npg("nrc", alpha = 0.7)(9)
######################### 20240830
batch_f<- list(sciHiC_gm12878="zhichao_biotinGM12878_20241012",
               dip_C_gm12878="zhichao_gm12878dipC_20241015",
               snHiChew_gm12878="zhichao_snHiChew400cellgm12878_20240918")

batch_data <- map2(batch_f, names(batch_f), function(file, cell_lib) {
  dt <- fread(paste0(file, '_downsample_rmdup_validpair.xls'))  # Read the file into a data.table
  dt[, cell_lib := cell_lib]  # Add the label column
  return(dt)})

# Row bind the tables
combined_table <- rbindlist(batch_data)

colnames(combined_table)<-c('readNum_per_cell','mean_rmdup_valid_pair','sd','sample')
#stack_df<-stack_df[stack_df$readNum_per_cell<6400000,]
combined_table<-combined_table[combined_table$readNum_per_cell<8000000,]
p1<-ggplot(combined_table, aes(x = readNum_per_cell, y = mean_rmdup_valid_pair,group=sample,color=sample))+ 
  geom_line(alpha=0.8,lwd=1,aes(linetype=sample))+scale_linetype_manual(values=c(3,1,5,2))+ 
  geom_point()+  geom_pointrange( aes(x=readNum_per_cell, ymin=mean_rmdup_valid_pair-sd, ymax=mean_rmdup_valid_pair+sd), alpha=0.7, size=0.5)+
  theme_classic()+scale_color_manual(values=c(mypal[1],mypal[2],mypal[4],mypal[6]))+
  labs(title='rmd_validpair',subtitle='gradient: 25000 50000 100000 200000 400000 800000\nnagano2017 400 cells \ndipC 400 cells\nhichew 400 cells')




batch_f<- list(sciHiC_gm12878="zhichao_biotinGM12878_20241012",
               dip_C_gm12878="zhichao_gm12878dipC_20241015",
               snHiChew_gm12878="zhichao_snHiChew400cellgm12878_20240918")

batch_data <- map2(batch_f, names(batch_f), function(file, cell_lib) {
  dt <- fread(paste0(file, '_downsample_dup.xls'))  # Read the file into a data.table
  dt[, cell_lib := cell_lib]  # Add the label column
  return(dt)})

# Row bind the tables
combined_table <- rbindlist(batch_data)
colnames(combined_table)<-c('readNum_per_cell','mean_dup','sd','sample')
combined_table<-combined_table[combined_table$readNum_per_cell<8000000,]

p2<-ggplot(combined_table, aes(x = readNum_per_cell, y = mean_dup,group=sample,color=sample))+ 
  geom_line(alpha=0.8,lwd=1,aes(linetype=sample))+scale_linetype_manual(values=c(3,1,5,2))+ 
  geom_point()+  geom_pointrange( aes(x=readNum_per_cell, ymin=mean_dup-sd, ymax=mean_dup+sd), alpha=0.7, size=0.5)+
  theme_classic()+scale_color_manual(values=c(mypal[1],mypal[2],mypal[4],mypal[6]))+
  labs(title='validpair dup',subtitle='gradient: 25000 50000 100000 200000 400000 800000\nnagano2017 400 cells \ndipC 400 cells\nhichew 400 cells')
p2

library(gridExtra)
pdf('saturation_curve_validpair_dup_zhichao_gm12878.pdf',height=7,width=12)
grid.arrange(p1,p2,nrow=1,ncol=2)
dev.off()

svg('saturation_curve_validpair_dup_zhichao_gm12878.svg',height=7,width=12)
grid.arrange(p1,p2,nrow=1,ncol=2)
dev.off()
