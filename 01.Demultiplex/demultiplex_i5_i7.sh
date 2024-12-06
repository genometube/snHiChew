lane_id=V350286856_L01
lib_id=240923_M004_V350286856_L01_CWHPE24090158-N1
proj_id=zhichao_snHiChew400cellgm12878_20240930
proj_dir=/research/xieyeming1/proj_2023/zhichao_snHiChew400cellgm12878_20240930
cell_num=2000

raw_R1_fq_raw="${proj_dir}/raw_data/${lib_id}/${lane_id}_read_1.fq.gz" # PATH to raw R1 fastq.gz, must be in gz format
raw_R2_fq_raw="${proj_dir}/raw_data/${lib_id}/${lane_id}_read_2.fq.gz" # PATH to raw R2 fastq.gz, must be in gz format

sample_index1="sample_index1_20240913.csv"
sample_index2="sample_index2.csv"

top_cell_num="2000" # number of cells to demultiplex, at least few hundreds greater than wet-lab experiment estimates
sampling_read_num="1000000"  # initial subsampling of raw fastq to identify top XXX cell barcode to demultiplex
BUFFER_SIZE="10000000"  # lines of fastq to store in ram
THREAD="4"  # number of thread for multiprocessing
TOP_CELL_BARCODE_LIST="top_${top_cell_num}_cells.list"
python='/research/zhangchen/software/anaconda3/bin/python' # path to python3, tested on Python 3.8.6

################################################################################################################
echo "get top cells ranked by read numbers"
echo ==========START at `date`==========
${python} demultiplex_hichew_top_cell.py -m sample_index.meta_data -1 ${sample_index1} -2 ${sample_index2} \
-R1 ${raw_R1_fq_raw} -R2 ${raw_R2_fq_raw} -d ./ -s ${sampling_read_num} -t ${top_cell_num}


################################################################################################################
echo "demultiplex begin"
decho ==========START at `date`==========
mkdir -p chunk
time ${python} demultiplex_hichew_multicore_gz.py -m sample_index.meta_data -1 ${sample_index1} -2 ${sample_index2} \
-d chunk -R1 ${raw_R1_fq_raw} -R2 ${raw_R2_fq_raw} -b ${BUFFER_SIZE} -t ${THREAD} -l ${TOP_CELL_BARCODE_LIST}

echo ==========END at `date`==========
################################################################################################################
#  echo "optional: cp to target dir"|awk '{print "\n"$0}'
#  echo ==========START at `date`==========
#  target_dir="${proj_dir}/demultiplex/${lib_id}/top${top_cell_num}"
#
#  mkdir -p ${target_dir}
#  cp chunk/* ${target_dir}/
#  echo ==========END at `date`==========
################################################################################################################