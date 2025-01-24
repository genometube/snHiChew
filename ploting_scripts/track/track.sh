#
#  # batchfile=($(ls -d /research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/snHiChew_HEK293Twt/hic_pro/*/|rev|cut -c 2-|rev))
#  batchfile=($(cat /research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/snHiChew_HEK293Twt/benchmark/hic_pro_metrics_sum/sorted_snHiChew_HEK293Twt_cut.xls|awk '{print $1}'))
#
#  sample_id=(${batchfile[@]##*/})
#
#  sample="snHiChew_HEK293Twt"
#  echo ${sample_id[@]}|wc
#
#  date
#
#  hic_pro_dir="/research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/snHiChew_HEK293Twt/hic_pro"
#
#  for i in ${sample_id[@]}; do
#   cat ${hic_pro_dir}/${i}/outdir/hic_results/data/${i}/${i}.allValidPairs >> ${sample}.allValidPairs
#  done
#
#  hicpro2higlass="/research/xieyeming1/software/hic_pro/hic_pro_3.1.0/HiC-Pro_3.1.0/bin/utils/hicpro2higlass_xieyeming1.sh"
#  juicer_tools="/research/liwenqing/Software/juicer/CPU/common/juicer_tools.jar"
#  chrom_size="/research/xieyeming1/db/hic_pro_ref/hg19/chrom_hg19.sizes"
#  resfrag="/research/xieyeming1/db/hic_pro_ref/hg19/DpnII_resfrag_hg19.bed"
#  hicpro2juicebox="/research/xieyeming1/software/hic_pro/hic_pro_3.1.0/HiC-Pro_3.1.0/bin/utils/hicpro2juicebox.sh"
#  #  cooler="/research/zhangchen/software/anaconda3/envs/liwenqing/bin/cooler"
#  #  pairix="/research/xieyeming1/software/pairix/pairix-master/bin/pairix"
#  #  bgzip="/research/xieyeming1/software/pairix/pairix-master/bin/bgzip"
#
#  export PATH=$PATH:/research/zhangchen/software/anaconda3/envs/liwenqing/bin
#  # conda activate /research/zhangchen/software/anaconda3/envs/liwenqing/
#  source /research/zhangchen/software/anaconda3/bin/activate liwenqing
#
#  export PATH=$PATH:/research/xieyeming1/software/pairix/pairix-master/bin
#
#  mkdir -p tmp
#  bash ${hicpro2juicebox} -i ${sample}.allValidPairs -g ${chrom_size} -j ${juicer_tools} -r ${resfrag}
#  bash ${hicpro2higlass} -i ${sample}.allValidPairs -r 5000 -c ${chrom_size} -p 16 -n -t tmp

# rm snHiChew_HEK293Twt*_tad_score.bm snHiChew_HEK293Twt*_boundaries.gff snHiChew_HEK293Twt*_domains.bed snHiChew_HEK293Twt*_boundaries.bed snHiChew_HEK293Twt*_score.bedgraph

sample="snHiChew_HEK293Tctcfkd"

fdr="0.05"
cool="/research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/hichew_ctcfko/hichew_readID_validpair/combined_mcool_stack_bulk/hichew_ctcfko_chr7.cool"
hicFindTADs="/research/zhangchen/software/anaconda3/envs/liwenqing/bin/hicFindTADs"
# ${hicFindTADs} -m ${sample}.mcool::resolutions/5000 \
${hicFindTADs} -m ${cool} \
--outPrefix ${sample}_${fdr} \
--minDepth 75000 \
--maxDepth 400000 \
--step 75000 \
--thresholdComparisons ${fdr} \
--delta 0.01 \
--correctForMultipleTesting fdr \
-p 16

hicPlotTADs="/mnt/software/anaconda3/envs/hicexplorer/bin/hicPlotTADs"
${hicPlotTADs} --tracks tracks.ini_compare_wt_kd --region chr7:145000000-148000000  -o subTAD_calling_compare_wt_kd.png

# ${hicPlotTADs} --tracks tracks.ini --region chr7:5000000-10000000  -o TAD_calling_comparison.png
# ${hicPlotTADs} --tracks tracks.ini_compare_final --region chr7:145000000-155000000  -o subTAD_calling_1_fdr_compare_final.png

# ${hicPlotTADs} --tracks tracks.ini.test --region chr7:145000000-155000000  -o subTAD_calling_1.png
