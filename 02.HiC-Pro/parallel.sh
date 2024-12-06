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

# The parallel command below can be replaced by any flow control script or software to submit jobs in parallel.
echo ${sample_id[@]}|sed 's/\ /\n/g'|parallel "sh hic-pro-single.sh {}"

