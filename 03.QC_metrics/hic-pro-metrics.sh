# proj_dir=/research/xieyeming1/proj_2024/luomeng_deepSnHichew_20241127
echo ==========START at `date`==========
bedtools="/research/xieyeming1/software/Miniconda/envs/xieyeming1/bin/bedtools" # path to bedtools, tested on v2.31.1
gatc_bed="/research/xieyeming1/db/hic_pro_ref/hg19/hg19_gatc.bed" # not required
frag_bed="/research/xieyeming1/db/hic_pro_ref/hg19/DpnII_resfrag_hg19.bed" # typical hic pro reference files for HiChew data

mkdir -p out_bed out_xls
sample="$1"
hic_results="${proj_dir}/hic_pro/${lib_id}/${sample}/outdir/hic_results" # PATH to hic pro results for each cell

date

cat ${hic_results}/data/${sample}/${sample}_*.fa.bwt2pairs.DEPairs|awk '$2==$5'|cut -f2,3,6 > out_bed/${sample}_DE.bed
${bedtools} intersect -a out_bed/${sample}_DE.bed -b ${gatc_bed} -wao|awk '{print $0"\t"$6-$2"\t"$3-$6}' > out_bed/${sample}_DE_overlap.bed

gatc_non=($(cat out_bed/${sample}_DE_overlap.bed|awk '$4=="."'|wc -l))
gatc_mid=($(cat out_bed/${sample}_DE_overlap.bed|awk '$4!="."'|awk '$8>5&&$9>5'|wc -l))
gatc_side=($(cat out_bed/${sample}_DE_overlap.bed|awk '$4!="."'|awk '$8<=5||$9<=5'|wc -l))

rm out_bed/${sample}_DE.bed

# genome fragment coverage
validpair_frag_num=($(cat ${hic_results}/data/${sample}/${sample}.allValidPairs|awk '{print $9"\n"$10}'|awk '{A[$1]++}END{for(i in A)print i,A[i]}' |wc -l))
all_frag_num=($(cat ${frag_bed}|cut -f4|awk '{A[$1]++}END{for(i in A)print i,A[i]}' |wc -l))
validpair_genome_frag_cov=($(awk "BEGIN {printf \"%.4f\n\", $validpair_frag_num / $all_frag_num }"))

# mapping
Total_pairs_processed=($(cat ${hic_results}/stats/${sample}/${sample}.mpairstat|grep -v "##"|sed -n '1p'|cut -f2))
Unmapped_pairs=($(cat ${hic_results}/stats/${sample}/${sample}.mpairstat|grep -v "##"|sed -n '2p'|cut -f2))
Low_qual_pairs=($(cat ${hic_results}/stats/${sample}/${sample}.mpairstat|grep -v "##"|sed -n '3p'|cut -f2))
Pairs_with_singleton=($(cat ${hic_results}/stats/${sample}/${sample}.mpairstat|grep -v "##"|sed -n '6p'|cut -f2))
Reported_pairs=($(cat ${hic_results}/stats/${sample}/${sample}.mpairstat|grep -v "##"|sed -n '10p'|cut -f2))
# pair filter
valid_pairs=($(cat ${hic_results}/stats/${sample}/${sample}.mRSstat|tail -n 11|head -1|cut -f2))
Dangling_end_pairs=($(cat ${hic_results}/stats/${sample}/${sample}.mRSstat|tail -n 6|head -1|cut -f2))
Religation_pairs=($(cat ${hic_results}/stats/${sample}/${sample}.mRSstat|tail -n 5|head -1|cut -f2))
Self_Cycle_pairs=($(cat ${hic_results}/stats/${sample}/${sample}.mRSstat|tail -n 4|head -1|cut -f2))
Dumped_pairs=($(cat ${hic_results}/stats/${sample}/${sample}.mRSstat|tail -n 1|head -1|cut -f2))
# dup, range
valid_interaction=($(cat ${hic_results}/stats/${sample}/${sample}_allValidPairs.mergestat|grep -v "##"|sed -n '1p'|cut -f2))
valid_interaction_rmdup=($(cat ${hic_results}/stats/${sample}/${sample}_allValidPairs.mergestat|grep -v "##"|sed -n '2p'|cut -f2))
trans_interaction=($(cat ${hic_results}/stats/${sample}/${sample}_allValidPairs.mergestat|grep -v "##"|sed -n '3p'|cut -f2))
cis_shortRange=($(cat ${hic_results}/stats/${sample}/${sample}_allValidPairs.mergestat|grep -v "##"|sed -n '5p'|cut -f2))
cis_longRange=($(cat ${hic_results}/stats/${sample}/${sample}_allValidPairs.mergestat|grep -v "##"|sed -n '6p'|cut -f2))

# mean_frag_size
mean_frag_size=($(cat ${hic_results}/data/${sample}/${sample}_*.fa.bwt2pairs.validPairs |cut -f8|awk '$1<100000'|awk '{ total += $1 } END { print total/NR }'))

echo "${sample} ${Total_pairs_processed} ${Unmapped_pairs} ${Low_qual_pairs} ${Pairs_with_singleton} ${Reported_pairs} ${Dangling_end_pairs} ${Religation_pairs} ${Self_Cycle_pairs} ${Dumped_pairs} ${valid_pairs} ${valid_interaction_rmdup} ${trans_interaction} ${cis_shortRange} ${cis_longRange} ${gatc_non} ${gatc_mid} ${gatc_side} ${validpair_genome_frag_cov} ${mean_frag_size}"| tr [:blank:] \\t > out_xls/${sample}_metrics.xls

# dup level sum
mkdir -p out_dupLevel
cat ${hic_results}/data/${sample}/${sample}_*.fa.bwt2pairs.validPairs |cut -f2-|sed 's/\t/\|/g'|awk '{A[$1]++}END{for(i in A)print i,A[i]}' |awk '{print $2}'|awk '{A[$1]++}END{for(i in A)print i,A[i]}'|awk '{print $0" "$1*$2}' |tr [:blank:] \\t| awk '{b[$1]=$3;sum=sum+$3} END{for (i in b) print i,b[i],(b[i]/sum)*100}'|tr [:blank:] \\t|sort -k1 -nr > out_dupLevel/${sample}_dup_freq_test.xls

echo ==========END at `date`==========
