# Please check image output and adjust knee1 knee2 and interval to get optimal cutoff based on the curve
# Example command: sh rm_empty_cell_auto.sh your_sample_id 1 2 150

sample="$1"
knee1_sensitivity="$2"
knee2_sensitivity="$3"
knee1_knee2_interval="$4"

# Example parameters:
# knee1_sensitivity=1
# knee2_sensitivity=2
# knee1_knee2_interval=150
# sample="zhichao_400cellgm12878b2_20240729"

Rscript='/mnt/software/anaconda3/envs/Bioinfo/bin/Rscript' # tested on R version 4.2.2, required R package can be found in rm_empty_cell_auto.R

${Rscript} rm_empty_cell_auto.R hic_pro_metrics.xls empty_filtered_${sample}_hic_pro_metrics.xls ${sample} \
 ${knee1_sensitivity} ${knee2_sensitivity} ${knee1_knee2_interval}

echo ==========END at `date`==========
