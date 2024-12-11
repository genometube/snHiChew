HiC_Pro="/research/xieyeming1/software/hic_pro/hic_pro_3.1.0/HiC-Pro_3.1.0/bin/HiC-Pro" # path to HiC-Pro
# lane_id=V350286856_L01
# lib_id=240923_M004_V350286856_L01_CWHPE24090158-N1
# proj_id=zhichao_snHiChew400cellgm12878_20240930
# cell_num=2000
# proj_dir=/research/xieyeming1/proj_2023/zhichao_snHiChew400cellgm12878_20240930

echo ==========START at `date`==========
fq_dir="${proj_dir}/demultiplex/${lib_id}/top${cell_num}" # PATH to demultiplexed fastq for each cell, must be in gz format
sample="$1"
mkdir -p ${sample}/rawdata/${sample}
ln -s ${fq_dir}/${sample}_R1.fq.gz ${proj_dir}/hic_pro/${lib_id}/${sample}/rawdata/${sample}/${sample}_R1.fastq
ln -s ${fq_dir}/${sample}_R2.fq.gz ${proj_dir}/hic_pro/${lib_id}/${sample}/rawdata/${sample}/${sample}_R2.fastq

cd ${sample}
${HiC_Pro} -i rawdata/ -o outdir -c ../config-hicpro.txt  # typical hic pro reference files for HiChew data

rm -rf outdir/bowtie_results/

echo ==========END at `date`==========
