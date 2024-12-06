Rscript="/mnt/software/anaconda3/envs/Bioinfo/bin/Rscript"
ins=/research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/snHiChew_HEK293Twt/EDA/validpair_downsample_mcool/validpair_500000/hichew_293t_500000_ins_cell_sum_25000_hg19.xls

length=($(cat hek293_gene_exp.xls|sed 's/.$//'|awk '{print $0"\t"$3-$2}'|awk '$6>300000'|wc -l))
long_gene=($(cat hek293_gene_exp.xls|sed 's/.$//'|awk '{print $0"\t"$3-$2}'|awk '$6>300000'|cut -f1-5))

echo ${long_gene[@]}
echo ${long_gene[@]}|wc

for i in $(seq 1 $length); do
  echo ${i}
  mkdir -p melt_score/gene_${i}
  tmp=($(awk "BEGIN {printf \"%.0f\n\", $i * 5 }"))
  chr=${long_gene[$tmp-5]}
  start=${long_gene[$tmp-4]}
  end=${long_gene[$tmp-3]}
  exp=${long_gene[$tmp-2]}
  gene=${long_gene[$tmp-1]}

  cat $ins|awk -v a="${chr}" '$1==a'|awk -v a="${start}" '$2>a'|awk -v a="${end}" '$3<a'|awk -v a="${gene}" '{print $0"\t"a}' > melt_score/gene_${i}/gene_${i}_ins.xls

  ${Rscript} melting_score.R melt_score/gene_${i}/gene_${i}_ins.xls melt_score/gene_${i}/gene_${i}_melt_score.xls \
  melt_score/gene_${i}/gene_${i}_cell_cluster_1.xls melt_score/gene_${i}/gene_${i}_cell_cluster_2.xls

done
