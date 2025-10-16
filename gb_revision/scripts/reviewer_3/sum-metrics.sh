date
cat out_xls/*xls > hic_pro_metrics.xls.tmp
echo "sample Total_pairs_processed Unmapped_pairs Low_qual_pairs Pairs_with_singleton Reported_pairs Dangling_end_pairs Religation_pairs Self_Cycle_pairs Dumped_pairs valid_pairs valid_interaction_rmdup trans_interaction cis_shortRange cis_longRange gatc_non gatc_mid gatc_side validpair_genome_frag_cov mean_frag_size"| tr [:blank:] \\t > hic_pro_metrics.header
cat hic_pro_metrics.header hic_pro_metrics.xls.tmp > hic_pro_metrics.xls

Rscript="/research/xieyeming1/software/Miniconda/envs/fyt/bin/Rscript"
${Rscript} hic_pro_metrics.R

date