#  batchfile=($(find ${fq_dir} -name "*_suffix" -type f))
#  sample_id_tmp=(${batchfile[@]##*/})
#  sample_id=(${sample_id_tmp[@]%_suffix})
#  echo ${sample_id[@]}
#
#  echo ${sample_id[@]}|sed 's/\ /\n/g'|parallel "sh melti-score.sh {}"
zcat /research/xieyeming1/db/attribute_table/hg19_rna_seq.bed.gz|awk '{print $0"\t"$3-$2}'|awk '$9>300000'|cut -f1-4,6,7 > hek293_long_tx_exp.xls

length=($(cat hek293_long_tx_exp.xls|wc -l))
echo $(seq 1 $length)|sed 's/\ /\n/g'|parallel "sh melti-score.sh {}"


