lane_id=V350278413_L03
lib_id=240726_M006_V350278413_L03_CWHPE24071633-N5
proj_dir=/research/xieyeming1/proj_2023/zhichao_400cellgm12878b1_20240729
proj_id=zhichao_400cellgm12878b1_20240729
top_cell_num=2000
#  proj_id="luomeng_deepSnHichew_20240428"
#  proj_path="/research/xieyeming1/proj_2024/${proj_id}"

# lib_id_raw=($(ls -d ${proj_dir}/raw_data/*M*_V*_L*))
# lib_id=(${lib_id_raw[@]##*/})
# lane_id=($(echo ${lib_id[@]}|cut -d "_" -f3,4))
raw_R1_fq_raw="${proj_dir}/raw_data/${lib_id}/${lane_id}_read_1.fq.gz"
raw_R2_fq_raw="${proj_dir}/raw_data/${lib_id}/${lane_id}_read_2.fq.gz"

sample_index1="sample_index1.csv"
sample_index2="sample_index2.csv"


sampling_read_num="10000000"
BUFFER_SIZE="20000000"
THREAD="8"
TOP_CELL_BARCODE_LIST="top_${top_cell_num}_cells.list"

python='/research/zhangchen/software/anaconda3/bin/python'
################################################################################################################
${python} demultiplex_hichew_top_cell.py -m sample_index.meta_data -1 ${sample_index1} -2 ${sample_index2} -R1 ${raw_R1_fq_raw} -R2 ${raw_R2_fq_raw} -d ./ -s ${sampling_read_num} -t ${top_cell_num}

echo "demultiplex begin"
date

mkdir -p chunk
time ${python} demultiplex_hichew_multicore_gz.py -m sample_index.meta_data -1 ${sample_index1} -2 ${sample_index2} -d chunk -R1 ${raw_R1_fq_raw} -R2 ${raw_R2_fq_raw} -b ${BUFFER_SIZE} -t ${THREAD} -l ${TOP_CELL_BARCODE_LIST}

date

target_dir="${proj_dir}/demultiplex/${lib_id}/top${top_cell_num}"

mkdir -p ${target_dir}
date

################################################################################################################
echo "cp to target dir"
# Target directory
cp chunk/* ${target_dir}/

date

################################################################################################################
echo "upload"

ossutil64="/research/xieyeming1/archive/ossutil/ossutil64"
config_file="/research/xieyeming1/archive/ossutil/bgionline.cfg"
bgionline_dir="oss://bgionline-chenjing1023-new/xieyeming1"

date

${ossutil64} cp -r -f --config-file ${config_file} ${target_dir} ${bgionline_dir}/${proj_id}/demultiplex/${lib_id}/

date
