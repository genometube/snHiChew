lane_id=V350278413_L03
lib_id=240726_M006_V350278413_L03_CWHPE24071633-N5
proj_id=zhichao_400cellgm12878b1_20240729
cell_num=2000
proj_dir=/research/xieyeming1/proj_2023/zhichao_400cellgm12878b1_20240729
echo ==========START at `date`==========
#  proj_id="zhichao_porecHiC_20240318_ctcfko"
#  proj_dir="/research/xieyeming1/proj_2024/${proj_id}"
#  lane_id="V350228045_L04"
#  lib_id="240403_M024_${lane_id}_CWHPE24030598-N1"
#  cell_num="2000"

fq_dir="${proj_dir}/demultiplex/${lib_id}/top${cell_num}"
sample="$1"
mkdir -p ${sample}/rawdata/${sample}
ln -s ${fq_dir}/${sample}_R1.fq.gz ${proj_dir}/hic_pro/${lib_id}/${sample}/rawdata/${sample}/${sample}_R1.fastq
ln -s ${fq_dir}/${sample}_R2.fq.gz ${proj_dir}/hic_pro/${lib_id}/${sample}/rawdata/${sample}/${sample}_R2.fastq

cd ${sample}
/research/xieyeming1/software/hic_pro/hic_pro_3.1.0/HiC-Pro_3.1.0/bin/HiC-Pro -i rawdata/ -o outdir -c ../config-hicpro.txt
rm -rf outdir/bowtie_results/

echo ==========END at `date`==========
