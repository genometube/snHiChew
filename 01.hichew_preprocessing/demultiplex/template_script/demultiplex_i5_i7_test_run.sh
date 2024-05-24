#  lane_id="V350228045_L04"
#  lib_id="240403_M024_${lane_id}_CWHPE24030598-N1"
#  proj_id="snHiChew_20240318"
#  cell_num="2000"

raw_R1_fq_raw="/research/xieyeming1/proj_2023/${proj_id}/raw_data/${lib_id}/${lane_id}_read_1.fq.gz"
raw_R2_fq_raw="/research/xieyeming1/proj_2023/${proj_id}/raw_data/${lib_id}/${lane_id}_read_2.fq.gz"
sample_index1="sample_index1.csv"
sample_index2="sample_index2.csv"
demultiplex_10x_i1_i2_py="demultiplex_i1_i2_V3_single_cell.py"
python='/research/zhangchen/software/anaconda3/bin/python'
target_dir="/research/xieyeming1/proj_2023/${proj_id}/demultiplex/${lib_id}/top${cell_num}"
tmp_dir="./${lib_id}"

mkdir -p ${tmp_dir}
date
echo "combine R1 R2"
zcat ${raw_R1_fq_raw}|head -40000 > ${tmp_dir}/${lane_id}_R1.fq
zcat ${raw_R2_fq_raw}|head -40000 > ${tmp_dir}/${lane_id}_R2.fq
paste -d $'' ${tmp_dir}/${lane_id}_R1.fq ${tmp_dir}/${lane_id}_R2.fq > ${tmp_dir}/${lane_id}.combine.fq
rm ${tmp_dir}/${lane_id}_R*fq
date
mkdir -p chunk

echo "demultiplexing"
${python} ${demultiplex_10x_i1_i2_py} -m sample_index.meta_data -i ${tmp_dir}/${lane_id}.combine.fq -1 ${sample_index1} -2 ${sample_index2} -d chunk -l ${lib_id}
rm ${tmp_dir}/${lane_id}.combine.fq
date
