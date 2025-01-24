computeMatrix=/mnt/software/anaconda3/envs/fyt/bin/computeMatrix
plotHeatmap=/mnt/software/anaconda3/envs/fyt/bin/plotHeatmap
# sample="CTCFchip_IP_hek_WT_replicate1.dedup.chr7"
sample="$1"

#  $computeMatrix reference-point \
#   -S ${sample}.bw \
#   -R /research/xieyeming1/db/attribute_table/ctcf.chr7.bed \
#   --referencePoint center \
#   -a 2000 -b 2000 --missingDataAsZero\
#   -out ${sample}.tab.gz

$plotHeatmap \
 -m ${sample}.tab.gz\
 -out ${sample}.pdf \
 --heatmapHeight 15  \
 --refPointLabel ctcf_site.center \
 --regionsLabel ctcf_site \
 --plotTitle ${sample} --zMin 0 --zMax 80 --colorMap Reds --yMax 80 --yMin 0

