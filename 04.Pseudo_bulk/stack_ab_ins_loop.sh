hicpro2higlass="/research/xieyeming1/software/hic_pro/hic_pro_3.1.0/HiC-Pro_3.1.0/bin/utils/hicpro2higlass_xieyeming1.sh"
juicer_tools="/research/liwenqing/Software/juicer/CPU/common/juicer_tools.jar"
chrom_size="/research/xieyeming1/db/hic_pro_ref/hg19/chrom_hg19.sizes"
resfrag="/research/xieyeming1/db/hic_pro_ref/hg19/DpnII_resfrag_hg19.bed"
hicpro2juicebox="/research/xieyeming1/software/hic_pro/hic_pro_3.1.0/HiC-Pro_3.1.0/bin/utils/hicpro2juicebox.sh"
export PATH=$PATH:/research/zhangchen/software/anaconda3/envs/liwenqing/bin
source /research/zhangchen/software/anaconda3/bin/activate liwenqing
export PATH=$PATH:/research/xieyeming1/software/pairix/pairix-master/bin
cooltools="/research/zhangchen/software/anaconda3/envs/liwenqing/bin/cooltools"

sample="snHiChew_400cellgm12878_20240930"
hic_pro_metrics="/research/xieyeming1/proj_2023/zhichao_snHiChew400cellgm12878_20240930/hic_pro_metrics_sum/240923_M004_V350286856_L01_CWHPE24090158-N1/hic_pro_metrics.xls"
compartment_bin=500000
gc_bin="/research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/hichew/compartment_tad/hg19.${compartment_bin}_bin.gc"
hic_pro_dir="/research/xieyeming1/proj_2023/zhichao_snHiChew400cellgm12878_20240930/hic_pro/240923_M004_V350286856_L01_CWHPE24090158-N1"

date
batchfile=($(cat ${hic_pro_metrics}|awk '$2>5000'|grep -v sample|awk '{print $1}'))
sample_id=(${batchfile[@]##*/})
echo ${sample_id[@]}|wc

for i in ${sample_id[@]}; do
 cat ${hic_pro_dir}/${i}/outdir/hic_results/data/${i}/${i}.allValidPairs >> ${sample}.allValidPairs
done

mkdir -p tmp
bash ${hicpro2juicebox} -i ${sample}.allValidPairs -g ${chrom_size} -j ${juicer_tools} -r ${resfrag}
bash ${hicpro2higlass} -i ${sample}.allValidPairs -r 5000 -c ${chrom_size} -p 16 -n -t tmp

${cooltools} call-compartments --reference-track ${gc_bin} --n-eigs 1 ${sample}.mcool::resolutions/${compartment_bin} -o ${sample}.mcool.res_${compartment_bin}.eig
${cooltools} diamond-insulation -o ${sample}.mcool.res_50k.window_150k.ins ${sample}.mcool::resolutions/50000 150000

date
