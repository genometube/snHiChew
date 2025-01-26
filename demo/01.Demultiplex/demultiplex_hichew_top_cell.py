import argparse
from operator import itemgetter
import os
import gzip

parser = argparse.ArgumentParser(description='')
parser.add_argument('-R1', '--R1', required=True, help = '')
parser.add_argument('-R2', '--R2', required=True, help = '')
parser.add_argument('-1', '--sample_index1', required=True, help = '')
parser.add_argument('-2', '--sample_index2', required=True, help = '')
parser.add_argument('-m', '--meta_data', required=True, help = '')
parser.add_argument('-d', '--out_dir', required=True, help = '')
parser.add_argument('-s', '--sampling_read_number', required=True, help = '')
parser.add_argument('-t', '--top_cell_threshold', required=True, help = 'top XXX cells sorted by demultiplexed read number')

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
        i7_indexNum_dict[tmp[i]] = tmp[0]
        i7index_set.add(tmp[i])
read_id_w_tag_set= set()

barcode_counter_dict = {}
line_count = 0

with gzip.open(args.R1,'rt') as f1,gzip.open(args.R2,'rt') as f2:
    for line1,line2 in zip(f1,f2):
        if line_count >= int(args.sampling_read_number) * 4:
            break
        if line_count % 4 == 0:
            read_id = line1.strip()[1:29]
        if line_count % 4 == 1:
            seq = line1.strip()+line2.strip()
            sample_index1, i7index = seq[200:206], seq[210:218]
        if line_count % 4 == 3:
            qual = line1.strip()+line2.strip()
            R1_insert_qual, R2_insert_qual, R1_insert, R2_insert = qual[3:100], qual[103:200], seq[3:100], seq[103:200]
            if sample_index1 in sample_index1_set and i7index in i7index_set:
                if sample_index1 + '_' + i7index in barcode_counter_dict:
                    barcode_counter_dict[sample_index1 + '_' + i7index] += 1
                else:
                    barcode_counter_dict[sample_index1 + '_' + i7index] = 1

                demultiplex_count[meta_dict[index1_sample_dict[sample_index1],i7_indexNum_dict[i7index]]] += 1
        line_count += 1

# Sort the dictionary by values (counts) in descending order and retrieve the top 2000
top_sequences = sorted(barcode_counter_dict.items(), key=lambda item: item[1], reverse=True)[:int(args.top_cell_threshold)]

# Print the sequences along with their counts
with open(os.path.abspath(args.out_dir)+ '/top_' + args.top_cell_threshold + '_cells.list','w') as f:
    for sequence, count in top_sequences:
        f.write('\t'.join([sequence,str(count)]) + '\n')

total_demultiplex_count = 0
for k in demultiplex_count:
    print("sample_id: "+k)
    print(demultiplex_count[k])
    total_demultiplex_count += demultiplex_count[k]

print("total demultiplex read count")
print(total_demultiplex_count)
print("total read count")
print(int(line_count/4))

