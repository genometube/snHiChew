Rscript="/mnt/software/anaconda3/envs/Bioinfo/bin/Rscript"
ins=/research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/snHiChew_HEK293Twt/EDA/validpair_downsample_mcool/validpair_500000/hichew_293t_500000_ins_cell_sum_25000_hg19.xls

i="$1"
echo ${i}
mkdir -p melt_score/gene_${i}
long_gene=($(awk -v i="$i" 'NR == i' hek293_long_tx_exp.xls))
chr=${long_gene[0]}
start=${long_gene[1]}
end=${long_gene[2]}
exp=${long_gene[4]}
tx=${long_gene[3]}
gene=${long_gene[5]}

cat $ins|awk -v a="${chr}" '$1==a'|awk -v a="${start}" '$2>a'|awk -v a="${end}" '$3<a'|awk -v a="${tx}" '{print $0"\t"a}' > melt_score/gene_${i}/gene_${i}_ins.xls

${Rscript} melting_score.R melt_score/gene_${i}/gene_${i}_ins.xls melt_score/gene_${i}/gene_${i}_melt_score.xls \
melt_score/gene_${i}/gene_${i}_cell_cluster_1.xls melt_score/gene_${i}/gene_${i}_cell_cluster_2.xls


