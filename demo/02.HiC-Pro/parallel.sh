lane_id=V350278413_L03
lib_id=240726_M006_V350278413_L03_CWHPE24071633-N5
proj_id=zhichao_400cellgm12878b1_20240729
cell_num=2000
proj_dir=/research/xieyeming1/proj_2023/zhichao_400cellgm12878b1_20240729

fq_dir="${proj_dir}/demultiplex/${lib_id}/top${cell_num}" # PATH to demultiplexed single cell fastq, must be in gz format
batchfile=($(find ${fq_dir} -name "*_R1.fq.gz" -type f))
sample_id_tmp=(${batchfile[@]##*/})
sample_id=(${sample_id_tmp[@]%_R1.fq.gz})
echo ${sample_id[@]}

# The parallel command below can be replaced by any flow control script or software to submit jobs in parallel.
echo ${sample_id[@]}|sed 's/\ /\n/g'|parallel "sh hic-pro-single.sh {}"

