bedGraphToBigWig='/research/zhangchen/software/bedGraphToBigWig'
cat /research/xieyeming1/db/genome/hg19/Sequence/bwa_index_hg19mt/genome.fa.fai|cut -f1-2 > hg19.fai
ref_fai='hg19.fai'

eigen_ins_dir="/research/xieyeming1/proj_2023/hichew_paper_20230710/scHiC_benchmark/hichew_ctcfko/hichew_readID_validpair/combined_mcool_stack_bulk"
sample="snHiChew_HEK293Tctcfkd"
# cat ${eigen_ins_dir}/snHiChew_HEK293Twt.mcool.res_500000.eig.cis.vecs.tsv|cut -f1,2,3,6|grep -v chrom|awk 'BEGIN{FS=OFS="\t"} {if($4=="") $4=0; print}'|sort -k1,1 -k2,2n > ${sample}_eigen.bedgraph
# ${bedGraphToBigWig} ${sample}_eigen.bedgraph ${ref_fai} ${sample}_eigen.bw

cat ${eigen_ins_dir}/hichew_ctcfko_chr7.mcool.res_25000.window_100k.ins|cut -f1,2,3,5|grep -v chrom|awk 'BEGIN{FS=OFS="\t"} {if($4=="nan") $4=0; print}'|sort -k1,1 -k2,2n > ${sample}_ins.bedgraph
${bedGraphToBigWig} ${sample}_ins.bedgraph ${ref_fai} ${sample}_ins.bw

