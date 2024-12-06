# qsub -cwd -l vf=4G,num_proc=4 -P RNAProj -binding linear:4 -q bc_rd.q,bc.q -e log -o log [bash.sh]
date

sample="snHiChew_HEK293Twt"
mkdir -p tmp
# rmdvalidPairs_dir="/research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/snHiChew_testis/hic_pro/${batch}_fq"
# batchfile=($(ls -d /research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/${sample}/hic_pro/*/|rev|cut -c 2-|rev))
batchfile=($(cat /research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/snHiChew_HEK293Twt/benchmark/hic_pro_metrics_sum/filtered_snHiChew_HEK293Twt_hic_pro_metrics.xls|grep -v sample|cut -f1))
sample_id=(${batchfile[@]##*/})
hic_pro_dir="/research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/${sample}/hic_pro"
echo ${sample_id[@]}|wc
echo "cell_name cell_id chrom1 pos1 chrom2 pos2 count"|tr [:blank:] \\t > data.txt
cell_id=0
for i in ${sample_id[@]}; do
  echo ${i}
  cat ${hic_pro_dir}/${i}/outdir/hic_results/data/${i}/${i}.allValidPairs|awk '$2==$5'|\
  awk -v a="${i}A" -v b="${cell_id}" '{print a"_"b"\t"$2"\t"$3"\t"$5"\t"$6}'|\
  sed 's/\_/\t/g'|sed 's/\t/\_/'|awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t1"}' >> data.txt
  # echo "${i} ${cell_id} testis"|tr [:blank:] \\t >> label_info.txt
  ((cell_id++))
done

#  batch="testis_replicate_2"
#  # rmdvalidPairs_dir="/research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/snHiChew_testis/hic_pro/${batch}_fq"
#  sample_id=($(cat /research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/snHiChew_testis/benchmark/hic_pro_metrics_sum/sorted_${batch}_cut.xls|grep -v sample|cut -f1))
#  hic_pro_dir="/research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/${sample}/hic_pro"
#
#  for i in ${sample_id[@]}; do
#    echo ${i}
#    cat ${hic_pro_dir}/${batch}_fq/${i}/outdir/hic_results/data/${i}/${i}.allValidPairs|awk '$2==$5'|\
#    awk -v a="${i}B" -v b="${cell_id}" '{print a"_"b"\t"$2"\t"$3"\t"$5"\t"$6}'|\
#    sed 's/\_/\t/g'|sed 's/\t/\_/'|awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t1"}' >> data.txt
#    # echo "${i} ${cell_id} testis"|tr [:blank:] \\t >> label_info.txt
#    ((cell_id++))
#  done

# gzip -c data.txt > data.txt.gz
python="/research/xieyeming1/software/Miniconda/envs/xieyeming1/bin/python"
${python} label_info_pickle.py -i repli_score/label_info.txt -o label_info.pickle
