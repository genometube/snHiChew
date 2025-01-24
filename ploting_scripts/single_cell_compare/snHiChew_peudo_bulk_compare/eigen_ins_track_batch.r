library(ggplot2)
library(plyr)
library(reshape2)
library(zoo)
library(pheatmap)
#library(HiCcompare)
library(Matrix)
library(data.table)
library(ComplexHeatmap)

library(ggplot2)
library(purrr)
library(cowplot)

################################################################################################ a/b compartment track
# Create example data
batch_f<- list(biotin_bulk_HEK293T="hek293_ngsHiC_700M",hek293t200cellb1="snHiChew_HEK293Twt",hek293t200cellb2="snHiChew_200cell293T",
            gm12878_400cellb1="snHiChew_400cellgm12878_20240918",gm12878_400cellb2="snHiChew_400cellgm12878_20240930",bulk_gm12878_4dn="4dn_bulk_gm12878")
#batch_data<- map(batch_f,~fread(paste0(.x,'.mcool.res_500000.eig.cis.vecs.tsv')))
batch_data <- map(batch_f, ~ {
  data <- fread(paste0(.x, '.mcool.res_500000.eig.cis.vecs.tsv'))
  data$E1[is.na(data$E1)] <- 0
  return(data)
})
colors <- c("darkred", "darkred", "darkred","turquoise4", "turquoise4", "turquoise4")

# Function to create a line plot for a given batch of data
create_line_plot <- function(batch_data, batch_name,color) {
  df_chr7<-batch_data[batch_data$chrom=="chr7"&batch_data$end<=40000000&batch_data$start>=10000000,]
  ggplot(df_chr7, aes(x = as.integer(rownames(df_chr7)), y = E1)) + 
    geom_area(fill = color,alpha=0.9) + 
    theme_classic() + 
    coord_cartesian(expand = FALSE) + ggtitle(batch_name)+
    labs(x = "Genomic position", y = "Eigenvector value")
}

# Use purrr::map to create a list of plots
plots <- map2(batch_data, names(batch_data), ~create_line_plot(.x, .y, colors[which(names(batch_data) == .y)]))

# Combine the plots using cowplot::plot_grid
combined_plot <- plot_grid(plotlist = plots, ncol = 1)

# add title
title_plot <- ggdraw() + draw_label("a/b compartment value, chr7:10-40M, 500k bin", fontface = 'bold', x = 0.5, hjust = 0.5)

# Combine the title and the combined plot
final_plot <- plot_grid(title_plot, combined_plot, ncol = 1, rel_heights = c(0.1, 1))

pdf('snHiChew_cor_eigen_track.pdf',height = 10,width=6)
print(final_plot)
dev.off()

################################################################################################ ins track
batch_f<- list(biotin_bulk_HEK293T="hek293_ngsHiC_700M",hek293t200cellb1="snHiChew_HEK293Twt",hek293t200cellb2="snHiChew_200cell293T",
               gm12878_400cellb1="snHiChew_400cellgm12878_20240918",gm12878_400cellb2="snHiChew_400cellgm12878_20240930",bulk_gm12878_4dn="4dn_bulk_gm12878")
# batch_data<- map(batch_f,~fread(paste0(.x,'.mcool.res_50k.window_100k.ins')))
batch_data <- map(batch_f, ~ {
  data <- fread(paste0(.x, '.mcool.res_50k.window_100k.ins'))
  data$log2_insulation_score_100000[is.na(data$log2_insulation_score_100000)] <- 0
  return(data)
})
colors <- c("darkred", "darkred", "darkred","turquoise4", "turquoise4", "turquoise4")
batch_data
create_line_plot <- function(batch_data, batch_name,color) {
  df_chr7<-batch_data[batch_data$chrom=="chr7"&batch_data$start>=145000000&batch_data$end<=155000000,]
  ggplot(df_chr7, aes(x = as.integer(rownames(df_chr7)), y = log2_insulation_score_100000)) + 
    geom_area(fill = color,alpha=0.9) + theme_classic() + ylim(-2,2)+
    coord_cartesian(expand = FALSE) + ggtitle(batch_name)+
    labs(x = "Genomic position", y = "ins value")
}

plots <- map2(batch_data, names(batch_data), ~create_line_plot(.x, .y, colors[which(names(batch_data) == .y)]))
combined_plot <- plot_grid(plotlist = plots, ncol = 1)
title_plot <- ggdraw() + draw_label("insulation score, chr7:145-155M, 50k bin", fontface = 'bold', x = 0.5, hjust = 0.5)
final_plot <- plot_grid(title_plot, combined_plot, ncol = 1, rel_heights = c(0.1, 1))

pdf('snHiChew_cor_ins_track.pdf',height = 10,width=6)
print(final_plot)
dev.off()

################################################################################################ ins track



