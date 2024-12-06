lane_id=V350286856_L01
lib_id=240923_M004_V350286856_L01_CWHPE24090158-N1
proj_id=zhichao_snHiChew400cellgm12878_20240930
cell_num=2000
proj_dir=/research/xieyeming1/proj_2023/zhichao_snHiChew400cellgm12878_20240930

fq_dir="${proj_dir}/demultiplex/${lib_id}/top${cell_num}" # PATH to demultiplexed single cell fastq, must be in gz format
batchfile=($(find ${fq_dir} -name "*_R1.fq.gz" -type f))
sample_id_tmp=(${batchfile[@]##*/})
sample_id=(${sample_id_tmp[@]%_R1.fq.gz})
echo ${sample_id[@]}

echo ==========START at `date`==========
# The parallel command below can be replaced by any flow control script or software to submit jobs in parallel.
echo ${sample_id[@]}|sed 's/\ /\n/g'|parallel "sh hic-pro-metrics.sh {}"

cat out_xls/*xls > hic_pro_metrics.xls.tmp  # if empty line occurs, rerun the HiC-Pro for this cell,
# or remove the empty lines by this command instead: cat out_xls/*xls |awk '$2>0' > hic_pro_metrics.xls.tmp

echo "sample Total_pairs_processed Unmapped_pairs Low_qual_pairs Pairs_with_singleton Reported_pairs Dangling_end_pairs Religation_pairs Self_Cycle_pairs Dumped_pairs valid_pairs valid_interaction_rmdup trans_interaction cis_shortRange cis_longRange gatc_non gatc_mid gatc_side validpair_genome_frag_cov mean_frag_size"| tr [:blank:] \\t > hic_pro_metrics.header
cat hic_pro_metrics.header hic_pro_metrics.xls.tmp > hic_pro_metrics.xls
echo ==========END at `date`==========
