hicPlotDistVsCounts="/mnt/software/anaconda3/envs/hicexplorer/bin/hicPlotDistVsCounts"

# mcool1=/research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/hichew/mcool_dot_hic/hek293_ngsHiC_700M.mcool
mcool2=/research/xieyeming1/proj_2023/zhichao_200cell293Tb1_20240729/benchmark/stack_ab_ins_loop/snHiChew_200cell293Tb1.mcool
mcool3=/research/xieyeming1/proj_2023/zhichao_200cell293Tb2_20240729/benchmark/stack_ab_ins_loop/snHiChew_200cell293Tb2.mcool
mcool4=/research/xieyeming1/proj_2023/zhichao_snHiChew400cellgm12878_20240918/stack_ab_ins_loop/snHiChew_400cellgm12878_20240918.mcool
mcool5=/research/xieyeming1/proj_2023/zhichao_snHiChew400cellgm12878_20240930/stack_ab_ins_loop/snHiChew_400cellgm12878_20240930.mcool
# mcool6=/research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/bulk_biotin_hic/GM12878/stack_ab_ins_loop/4dn_bulk_gm12878.mcool

res=50000
${hicPlotDistVsCounts} -m ${mcool2}::resolutions/${res} \
${mcool3}::resolutions/${res} ${mcool4}::resolutions/${res} \
${mcool5}::resolutions/${res} \
-o snHiChew_HEK_gm12878_distance_decay_${res}_no_bulk.pdf --labels 200cell293Tb1 200cell293Tb2 400cellgm12878b1 400cellgm12878b2 --maxdepth 300000000 --plotsize 5 4.2


#  mcool1=/research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/hichew/mcool_dot_hic/hek293_ngsHiC_700M.mcool
#  mcool2=/research/xieyeming1/proj_2023/zhichao_200cell293Tb1_20240729/benchmark/stack_ab_ins_loop/snHiChew_200cell293Tb1.mcool
#  mcool3=/research/xieyeming1/proj_2023/zhichao_200cell293Tb2_20240729/benchmark/stack_ab_ins_loop/snHiChew_200cell293Tb2.mcool
#  mcool4=/research/xieyeming1/proj_2023/zhichao_snHiChew400cellgm12878_20240918/stack_ab_ins_loop/snHiChew_400cellgm12878_20240918.mcool
#  mcool5=/research/xieyeming1/proj_2023/zhichao_snHiChew400cellgm12878_20240930/stack_ab_ins_loop/snHiChew_400cellgm12878_20240930.mcool
#  mcool6=/research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/bulk_biotin_hic/GM12878/stack_ab_ins_loop/4dn_bulk_gm12878.mcool
#
#  res=100000
#  ${hicPlotDistVsCounts} -m ${mcool1}::resolutions/${res} ${mcool2}::resolutions/${res} \
#  ${mcool3}::resolutions/${res} ${mcool4}::resolutions/${res} \
#  ${mcool5}::resolutions/${res} ${mcool6}::resolutions/${res} \
#  -o snHiChew_HEK_gm12878_distance_decay_${res}.pdf --labels biotin_bulk_HEK293T 200cell293Tb1 200cell293Tb2 400cellgm12878b1 400cellgm12878b2 4dn_bulk_gm12878  --maxdepth 300000000 --plotsize 5 4.2


# mcool1=/research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/snHiChew_HEK293Twt/benchmark/stack_ab_ins_loop/snHiChew_HEK293Twt.mcool
