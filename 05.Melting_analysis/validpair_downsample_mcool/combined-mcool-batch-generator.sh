thres=500000

thread="$1"

sample_id=($(cat /research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/snHiChew_HEK293Twt/benchmark/hic_pro_metrics_sum/filtered_snHiChew_HEK293Twt_hic_pro_metrics.xls|grep -v sample|awk -v t=${thres} '$12>t'|cut -f1))
echo ${sample_id[@]}|wc

total_num=${#sample_id[@]}
min=0
max=${#sample_id[@]}
max=($(awk "BEGIN {printf \"%.0f\n\", $total_num - 1 }"))
interval=($(awk "BEGIN {printf \"%.0f\n\", ($total_num / $thread)+1 }"))

echo "batch_num total_sample_num interval start_num end_num"|tr [:blank:] \\t > run_metrics.xls
echo "$thread $total_num $interval $min $max"|tr [:blank:] \\t >> run_metrics.xls
echo "batch_num:$thread total_sample_num:$total_num interval:$interval start_num:$min end_num:$max"|tr [:blank:] \\t

date

for i in $(seq 1 $thread); do
  batch_start=($(awk "BEGIN {printf \"%.0f\n\", $min + $interval * ($i - 1) }"))
  batch_end=($(awk "BEGIN {printf \"%.0f\n\", $min + $interval * $i - 1 }"))
  echo "sh interval_${i}.sh &" >> batch_run.sh
  for j in $(seq $batch_start $batch_end); do
    if [ $j -le $max ];then
      if test ! -d "${i}";then
        echo "sh combined-mcool.sh ${sample_id[j]}" >> interval_${i}.sh
      fi
    fi
  done
done

date