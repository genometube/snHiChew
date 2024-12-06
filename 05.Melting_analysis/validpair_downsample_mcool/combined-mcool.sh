thres=500000
cooltools="/research/zhangchen/software/anaconda3/envs/liwenqing/bin/cooltools"

proj_dir="/research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/snHiChew_HEK293Twt/EDA/validpair_downsample_mcool"

compartment_bin=1000000
genome_gc="/research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/hichew/hichew_readID_validpair/combined_mcool_chr7_validpair_downsample/hg19.${compartment_bin}_bin.gc"
thread=2

sample="$1"
wd="${proj_dir}/validpair_${thres}/mcool/${sample}"

mkdir -p ${wd}
cd ${wd}

chrom_size="/research/xieyeming1/db/hic_pro_ref/hg19/chrom_hg19.sizes"
hic_pro_dir="/research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/snHiChew_HEK293Twt/hic_pro"

date
hicpro2higlass="/research/xieyeming1/software/hic_pro/hic_pro_3.1.0/HiC-Pro_3.1.0/bin/utils/hicpro2higlass_xieyeming1.sh"

export PATH=$PATH:/research/zhangchen/software/anaconda3/envs/liwenqing/bin
export PATH=$PATH:/research/xieyeming1/software/pairix/pairix-master/bin
python=/research/xieyeming1/software/Miniconda/envs/xieyeming1/bin/python

cat ${hic_pro_dir}/${sample}/outdir/hic_results/data/${sample}/${sample}.allValidPairs|shuf|head -${thres} > ${sample}.allValidPairs
bash ${hicpro2higlass} -i ${sample}.allValidPairs -r 5000 -c ${chrom_size} -p ${thread} -n -t tmp

${cooltools} call-compartments --reference-track ${genome_gc} --n-eigs 1 ${sample}.mcool::resolutions/${compartment_bin} -o ${sample}.mcool.res_${compartment_bin}.eig

hicFindTADs="/research/zhangchen/software/anaconda3/envs/liwenqing/bin/hicFindTADs"
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
cd ${proj_dir}

date

#    mkdir -p tmp
#    cat /research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/hichew/hichew_readID_validpair/all_validpair/combined_pairs_20240405/combined_rmdvalidpair/${sample}.rmdvalidpairs|shuf|head -${thres}|sort -T tmp -k2,2 -k4,4 -k3,3n -k5,5n >  ${sample}.pairs
#    ${bgzip} ${sample}.pairs
#    ${pairix} -p pairs -f ${sample}.pairs.gz
#
#    # cat ${ref_fai}|awk '{print $1"\t"$2}'|grep -v chrM > chrom_size.txt
#
#    pairs_file=${sample}.pairs.gz  # bgzipped, with .px2 index
#    chrsize_file=${chr_size}
#    bin_size=50000
#    out_prefix=${sample}
#    max_split=2 # (default 2)
#
#    ${cooler} makebins -o $chrsize_file:$bin_size $chrsize_file $bin_size
#    ${cooler} cload pairix -p $ncores -s $max_split $chrsize_file:$bin_size $pairs_file $out_prefix.cool
#
#    chunksize=10000000
#    #RES_ARG="-r 25000,50000,100000,250000,500000,1000000,2500000"
#    RES_ARG="-r ${bin_size},${compartment_bin}"
#    ${cooler} zoomify --balance --balance-args '--convergence-policy store_nan' -n $ncores -o $out_prefix.mcool -c $chunksize $RES_ARG $out_prefix.cool
