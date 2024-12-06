date
python="/research/xieyeming1/software/Miniconda/envs/xieyeming1/bin/python"

#  validPairs_dir="/research/xieyeming1/proj_2023/hicVelocity_20231027/dowsample_293t/sc_read_validpair/validPairs_5000_bin_mat"
#  rmdvalidPairs_dir="${validPairs_dir}/rmdvalidPairs"
#  raw_dir="${validPairs_dir}/raw"
#  ice_dir="${validPairs_dir}/ice"
#  batchfile=($(ls ${rmdvalidPairs_dir}))
#  sample_id_tmp=(${batchfile[@]##*/})
#  sample_id=(${sample_id_tmp[@]%.rmdvalidPairs})
#  echo ${sample_id[@]}|wc

sample="snHiChew_HEK293Twt"
# batchfile=($(ls -d /research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/${sample}/hic_pro/*/|rev|cut -c 2-|rev))
batchfile=($(cat /research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/snHiChew_HEK293Twt/benchmark/hic_pro_metrics_sum/filtered_snHiChew_HEK293Twt_hic_pro_metrics.xls|grep -v sample|cut -f1))

sample_id=(${batchfile[@]##*/})

mkdir -p tmp data_tmp

# echo "cell_name cell_id chrom1 pos1 chrom2 pos2 count"|tr [:blank:] \\t > data.txt
cell_id=0
for i in ${sample_id[@]}; do
  echo ${i}
  # cat data_tmp/${i}_data.txt|awk -v a="${i}" -v b="${cell_id}" '{print a"\t"b"\t"$0}' >> data.txt
  cat tmp/${i}.early_ratio|awk -v a="${i}" -v b="${cell_id}" '{print a"\t"b"\t"$2}' >> label_info.txt
  ((cell_id++))
done

${python} label_info_pickle.py -i label_info.txt -o label_info.pickle

date
