# qsub -cwd -l vf=2G,num_proc=2 -P RNAProj -binding linear:2 -q bc_rd.q,bc.q -e log -o log [bash.sh]
# chrom="chr1"

date
compartment_bin=500000
cooltools="/research/zhangchen/software/anaconda3/envs/liwenqing/bin/cooltools"

${cooltools} genome binnify --all-names /research/xieyeming1/db/hic_pro_ref/hg19/chrom_hg19.sizes ${compartment_bin} > hg19.${compartment_bin}_bin
${cooltools} genome gc hg19.${compartment_bin}_bin /research/xieyeming1/db/genome/hg19/Sequence/bwa_index_hg19mt/genome.fa > hg19.${compartment_bin}_bin.gc

date
