date
#  lane_id="V350228045_L04"
#  lib_id="240403_M024_${lane_id}_CWHPE24030598-N1"
#  proj_id="snHiChew_20240318"
#  cell_num="2000"

thread="$1"
proj_dir="/research/xieyeming1/proj_2023/${proj_id}"

fq_dir="${proj_dir}/demultiplex/${lib_id}/top${cell_num}"
batchfile=($(find ${fq_dir} -name "*_R1.fq.gz" -type f))
sample_id_tmp=(${batchfile[@]##*/})
sample_id=(${sample_id_tmp[@]%_R1.fq.gz})
echo ${sample_id[@]}
total_num=${#sample_id[@]}

min=0
max=($(awk "BEGIN {printf \"%.0f\n\", $total_num  - 1 }"))
interval=($(awk "BEGIN {printf \"%.0f\n\", ($total_num / $thread)+1 }"))

echo "batch_num total_sample_num interval start_num end_num"|tr [:blank:] \\t > run_metrics.xls
echo "$thread $total_num $interval $min $max"|tr [:blank:] \\t >> run_metrics.xls
echo "batch_num:$thread total_sample_num:$total_num interval:$interval start_num:$min end_num:$max"|tr [:blank:] \\t


for i in $(seq 1 $thread); do
  batch_start=($(awk "BEGIN {printf \"%.0f\n\", $min + $interval * ($i - 1) }"))
  batch_end=($(awk "BEGIN {printf \"%.0f\n\", $min + $interval * $i - 1 }"))
  echo "sh interval_${i}.sh &" >> hic_pro_batch_run.sh
  for j in $(seq $batch_start $batch_end); do
    if [ $j -le $max ];then
      if test ! -d "${i}";then
        echo "sh hic-pro-single.sh ${sample_id[j]}" >> interval_${i}.sh
      fi
    fi
  done
done

date