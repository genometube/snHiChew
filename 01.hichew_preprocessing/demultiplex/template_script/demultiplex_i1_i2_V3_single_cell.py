import argparse
from operator import itemgetter
import os
# import gzip

parser = argparse.ArgumentParser(description='')
parser.add_argument('-i', '--input', required=True, help = '')
parser.add_argument('-1', '--sample_index1', required=True, help = '')
parser.add_argument('-2', '--sample_index2', required=True, help = '')
# parser.add_argument('-b', '--barcode_10x', required=True, help = '')
parser.add_argument('-m', '--meta_data', required=True, help = '')
parser.add_argument('-d', '--out_dir', required=True, help = '')
parser.add_argument('-l', '--lib', required=True, help = '')

args = parser.parse_args()

def rev_comp(DNAstr):
    cis = 'ATCGNatcgn'
    trans = 'TAGCNtagcn'
    rc = str.maketrans(cis,trans)
    rc_DNAstr = DNAstr.translate(rc)[::-1]
    return rc_DNAstr

with open(args.meta_data, 'r') as f:
    line = f.read().strip().split('\n')

meta_dict={}
demultiplex_count = {}
for i in range(len(line)):
    tmp = line[i].split(",")
    meta_dict[(tmp[1],tmp[2])] = tmp[0]
    demultiplex_count[tmp[0]] = 0

with open(args.sample_index1, 'r') as f:
    line = f.read().strip().split('\n')

sample_index1_set=set()
index1_sample_dict = {}
for i in range(len(line)):
    tmp = line[i].split(",")
    sample_index1_set.add(tmp[0])
    index1_sample_dict[tmp[0]] = tmp[1]
with open(args.sample_index2, 'r') as f:
    line = f.read().strip().split('\n')

i7_indexNum_dict = {}
i7index_set = set()
for i in range(len(line)):
    tmp = line[i].split(",")
    for i in range(1,len(tmp)):

        i7_indexNum_dict[(tmp[i])] = tmp[0]
        i7index_set.add((tmp[i]))

count = 0
with open(args.input) as f:
    for line in f:
        if count % 4 == 0:
            read_id = line.strip()[1:29]
        if count % 4 == 1:
            seq = line.strip()
            sample_index1, i7index = seq[200:206], seq[210:218]
        if count % 4 == 3:
            qual = line.strip()
            R1_insert_qual, R2_insert_qual, R1_insert, R2_insert = qual[3:100], qual[103:200], seq[3:100], seq[103:200]
            if i7index in i7index_set and sample_index1 in sample_index1_set:
                read_id_w_tag = '@' + read_id + '_' + sample_index1 + '_' + i7index

                with open(os.path.abspath(args.out_dir) + '/' + sample_index1 + '_' + i7index + '_R1.fq','a') as f:
                    f.write('\n'.join([read_id_w_tag,R1_insert,'+',R1_insert_qual])+'\n')

                with open(os.path.abspath(args.out_dir) + '/' + sample_index1 + '_' + i7index + '_R2.fq','a') as f:
                    f.write('\n'.join([read_id_w_tag, R2_insert, '+', R2_insert_qual]) + '\n')
                demultiplex_count[meta_dict[index1_sample_dict[sample_index1],i7_indexNum_dict[i7index]]] += 1
        count += 1

total_demultiplex_count = 0
for k in demultiplex_count:
    print(args.lib)
    # print(k)
    print(demultiplex_count[k])
    total_demultiplex_count += demultiplex_count[k]

print("total demultiplex count")
print(total_demultiplex_count)
