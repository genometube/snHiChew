date
# batchfile=($(ls -d /research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/snHiChew_HEK293Twt/hic_pro/*/|rev|cut -c 2-|rev))
batchfile=($(cat /research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/snHiChew_HEK293Twt/benchmark/hic_pro_metrics_sum/filtered_snHiChew_HEK293Twt_hic_pro_metrics.xls|grep -v sample|cut -f1))

sample_id=(${batchfile[@]##*/})
# sample_id=($(cat /research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/${sample}/benchmark/hic_pro_metrics_sum/hic_pro_metrics.xls|grep -v sample|cut -f1))

echo ${sample_id[@]}|wc

date

echo ${sample_id[@]}
for i in ${sample_id[@]}; do
  sh repli_score.sh ${i}
  # qsub -cwd -l vf=2G,num_proc=2 -binding linear:2 -P RNAProj -q bc_rd.q,bc.q -e log -o log repli_score.sh ${i}
done

date

#  cat mt_contact_near_ratio.xls.raw | \
#  awk '{if ($2<=0.5 && $3>=0.30)
#  print $0"\tPost_M";
#  else if ($2>0.5 && $2<=0.63)
#  print $0"\tG1";
#  else if ($2>0.63 && $2<=0.785)
#  print $0"\tEarly_to_mid_S";
#  else if ($2>0.785)
#  print $0"\tMid_S_to_G2";
#  else
#  print $0"\tPre_M"}' > mt_contact_near_ratio.xls




# echo ${sample_id[@]}|tr [:blank:] \\n|awk '{print $1"\t_"$1}' > sample.list
#  cat out_xls/*xls > hic_pro_metrics.xls.tmp
#  echo "sample Total_pairs_processed Unmapped_pairs Low_qual_pairs Pairs_with_singleton Reported_pairs Dangling_end_pairs Religation_pairs Self_Cycle_pairs Dumped_pairs valid_pairs valid_interaction_rmdup trans_interaction cis_shortRange cis_longRange Dangling_end_pairs_gatc_ratio validpair_genome_frag_cov mean_frag_size"| tr [:blank:] \\t > hic_pro_metrics.header
#  cat hic_pro_metrics.header hic_pro_metrics.xls.tmp > hic_pro_metrics.xls

