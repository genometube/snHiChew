date
# /research/xieyeming1/proj_2023/hichew_paper_20230710/repli_seq/hek293/E_L_ratio/E_L_293_replicate_1_20000Loess_overZero.bedGraph
bedtools="/research/xieyeming1/software/Miniconda/envs/xieyeming1/bin/bedtools"
batchfile=($(ls -d /research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/snHiChew_HEK293Twt/hic_pro/*/|rev|cut -c 2-|rev))
sample_id=(${batchfile[@]##*/})

# sample="snHiChew_HEK293Twt"
echo ${sample_id[@]}|wc

date

sample="$1"

mkdir -p tmp
hic_pro_dir="/research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/snHiChew_HEK293Twt/hic_pro"

cat ${hic_pro_dir}/${sample}/outdir/hic_results/data/${sample}/${sample}.allValidPairs|awk '{print $2"\t"$3-1"\t"$3"\n"$5"\t"$6-1"\t"$6}' > tmp/${sample}.bed

early_loci=($(${bedtools} intersect -a tmp/${sample}.bed -b E_L_293_replicate_1_20000Loess_overZero.bedGraph|wc -l))
all_loci=($(wc -l tmp/${sample}.bed))

early_ratio=($(awk "BEGIN {printf \"%.4f\n\", $early_loci / $all_loci }"))

echo "${sample} ${early_ratio}"|tr [:blank:] \\t > tmp/${sample}.early_ratio


# dist_binID=($(awk "BEGIN {printf \"%.0f\n\", (log(1116)/log(2)-10)/0.125 }"))

# 2^(10+0.125*143)
# (log2(1116)-10)/0.125
# bin=($(awk "BEGIN {printf \"%.0f\n\", (log(1116)/log(2)-10)/0.125 }"))

# dist_arr=($(cat ${hic_results}/data/${sample}/${sample}.allValidPairs | awk '$2==$5'|awk '{print $6-$3}'))
#  for i in ${dist_arr[@]}; do
#    awk "BEGIN {printf \"%.0f\n\", (log($i)/log(2)-10)/0.125 }"|awk '$1>37' >> tmp/${sample}.contact_dist
#  done

#  mt_contact=($(cat ${hic_results}/data/${sample}/${sample}.allValidPairs|awk '$2=="chrM"'|awk '$5=="chrM"'|wc -l))
#  all_contact=($(cat ${hic_results}/data/${sample}/${sample}.allValidPairs|wc -l))
#  mt_ratio=($(awk "BEGIN {printf \"%.2f\n\", $mt_contact / $all_contact }"))
