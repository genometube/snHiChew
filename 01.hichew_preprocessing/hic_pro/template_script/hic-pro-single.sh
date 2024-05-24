echo ==========START at `date`==========
#  lane_id="V350228045_L04"
#  lib_id="240403_M024_${lane_id}_CWHPE24030598-N1"
#  proj_id="snHiChew_20240318"
#  cell_num="2000"

proj_dir="/research/xieyeming1/proj_2023/${proj_id}"
fq_dir="${proj_dir}/demultiplex/${lib_id}/top${cell_num}"
sample="$1"

mkdir -p ${sample}/rawdata/${sample}
ln -s ${fq_dir}/${sample}_R1.fq.gz ${proj_dir}/hic_pro/${lib_id}/${sample}/rawdata/${sample}/${sample}_R1.fastq
ln -s ${fq_dir}/${sample}_R2.fq.gz ${proj_dir}/hic_pro/${lib_id}/${sample}/rawdata/${sample}/${sample}_R2.fastq

cd ${sample}
/research/xieyeming1/software/hic_pro/hic_pro_3.1.0/HiC-Pro_3.1.0/bin/HiC-Pro -i rawdata/ -o outdir -c ../config-hicpro.txt
rm -rf outdir/bowtie_results/

echo ==========END at `date`==========
