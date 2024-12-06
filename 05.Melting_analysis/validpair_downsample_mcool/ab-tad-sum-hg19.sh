thres=500000

# grep /research/xieyeming1/db/hic_pro_ref/hg19/chrom_hg19.sizes > hg19.sizes

bedtools="/research/xieyeming1/software/Miniconda/envs/xieyeming1/bin/bedtools"
compartment_bin=1000000
ins_bin=25000

$bedtools makewindows -g /research/xieyeming1/db/hic_pro_ref/hg19/chrom_hg19.sizes -w ${ins_bin} > hg19_${ins_bin}.bed

sample_id=($(cat /research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/snHiChew_HEK293Twt/benchmark/hic_pro_metrics_sum/filtered_snHiChew_HEK293Twt_hic_pro_metrics.xls|grep -v sample|awk -v t=${thres} '$12>t'|cut -f1))
echo ${sample_id[@]}|wc

sample="hichew_293t_${thres}"

date
for i in ${sample_id[@]}; do
  # cat mcool/${i}/${i}.mcool.res_${compartment_bin}.eig.cis.vecs.tsv|awk -v a="${i}" '{if ($5 != "") print i++"\t"$6"\t"a; else print i++"\t0\t"a}'|tail -n +2 >> ${sample}_eigen_cell_sum_${compartment_bin}.xls
  # cat mcool/${i}/${i}.mcool.res_50000.window_150k.ins|awk -v a="${i}" '{if ($5 != "nan") print i++"\t"$5"\t"a; else print i++"\t0\t"a}'|tail -n +2 >> ${sample}_ins_cell_sum.xls
  $bedtools intersect -a hg19_${ins_bin}.bed -b mcool/${i}/${i}_score.bedgraph -wao |awk '!seen[$1, $2, $3]++'|cut -f1-3,7|awk -v a="${i}" '{if ($4 == ".") print $1"\t"$2"\t"$3"\t0\t"a; else print $0"\t"a}' >> ${sample}_ins_cell_sum_${ins_bin}_hg19.xls
done

#  gzip -c ${sample}_eigen_cell_sum_${compartment_bin}.xls > ${sample}_eigen_cell_sum_${compartment_bin}.xls.gz
#  gzip -c ${sample}_ins_cell_sum_${ins_bin}.xls > ${sample}_ins_cell_sum_${ins_bin}.xls.gz

date
