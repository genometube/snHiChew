thres=500000 # valid pair thres

proj_dir="/research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/snHiChew_HEK293Twt/EDA/validpair_downsample_mcool"

compartment_bin=1000000
genome_gc="/research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/hichew/hichew_readID_validpair/combined_mcool_chr7_validpair_downsample/hg19.${compartment_bin}_bin.gc"
thread=2

chrom_size="/research/xieyeming1/db/hic_pro_ref/hg19/chrom_hg19.sizes" # hic pro GENOME_SIZE file
hic_pro_dir="/research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/snHiChew_HEK293Twt/hic_pro" # PATH to hic pro results for each cell

cooltools="/research/zhangchen/software/anaconda3/envs/liwenqing/bin/cooltools"
hicpro2higlass="/research/xieyeming1/software/hic_pro/hic_pro_3.1.0/HiC-Pro_3.1.0/bin/utils/hicpro2higlass_xieyeming1.sh"
export PATH=$PATH:/research/zhangchen/software/anaconda3/envs/liwenqing/bin
export PATH=$PATH:/research/xieyeming1/software/pairix/pairix-master/bin
python=/research/xieyeming1/software/Miniconda/envs/xieyeming1/bin/python
hicFindTADs="/research/zhangchen/software/anaconda3/envs/liwenqing/bin/hicFindTADs"

sample="$1" # cell barcode
wd="${proj_dir}/validpair_${thres}/mcool/${sample}"

mkdir -p ${wd}
cd ${wd}

date

cat ${hic_pro_dir}/${sample}/outdir/hic_results/data/${sample}/${sample}.allValidPairs|shuf|head -${thres} > ${sample}.allValidPairs
bash ${hicpro2higlass} -i ${sample}.allValidPairs -r 5000 -c ${chrom_size} -p ${thread} -n -t tmp

${cooltools} call-compartments --reference-track ${genome_gc} --n-eigs 1 ${sample}.mcool::resolutions/${compartment_bin} -o ${sample}.mcool.res_${compartment_bin}.eig

${hicFindTADs} -m ${sample}.mcool::resolutions/50000 \
--outPrefix ${sample} \
--minDepth 200000 \
--maxDepth 750000 \
--step 50000 \
--thresholdComparisons 0.05 \
--delta 0.01 \
--correctForMultipleTesting fdr \
-p ${thread}
#  rm ${sample}.cool ${sample}.allValidPairs chrom_size.txt*

date
