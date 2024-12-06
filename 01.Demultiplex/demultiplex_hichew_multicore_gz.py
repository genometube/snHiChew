import argparse
from operator import itemgetter
import os
import gzip
from multiprocessing import Process, Queue
import time

parser = argparse.ArgumentParser(description='')
# parser.add_argument('-i', '--input', required=True, help = '')
parser.add_argument('-R1', '--R1', required=True, help = 'R1 fastq.gz')
parser.add_argument('-R2', '--R2', required=True, help = 'R2 fastq.gz')
parser.add_argument('-1', '--sample_index1', required=True, help = '')
parser.add_argument('-2', '--sample_index2', required=True, help = '')
parser.add_argument('-m', '--meta_data', required=True, help = '')
parser.add_argument('-d', '--out_dir', required=True, help = '')
parser.add_argument('-b', '--buffer_size', required=True, help = '')
parser.add_argument('-t', '--thread', required=True, help = '')
parser.add_argument('-l', '--top_cell_barcode_list', required=True, help = '')

args = parser.parse_args()

def writeFile(pFileName, pReadArray):
    with gzip.open(pFileName, 'at') as f:
        for read in pReadArray:
            f.write(read)

def writeSampleFiles(pFileNameList, pBuffer, pQueue):
    for i, cells in enumerate(pFileNameList):
        for j, cell_ in enumerate(cells):
            writeFile(cell_, pBuffer[i][j])
    pQueue.put('Done')

def handleCompressingMulticore(pFileNameList, pBuffer, pThreads):
    filesPerThread = len(pFileNameList) // pThreads
    if filesPerThread == 0:
        writeSampleFiles(pFileNameList, pBuffer, None)
    else:
        queue = [None] * pThreads
        process = [None] * pThreads
        file_list_sample = [None] * pThreads
        buffer_sample = [None] * pThreads

        all_data_collected = False

        for i in range(pThreads):
            if i < pThreads - 1:
                file_list_sample = pFileNameList[i * filesPerThread:(i + 1) * filesPerThread]
                buffer_sample = pBuffer[i * filesPerThread:(i + 1) * filesPerThread]
            else:
                file_list_sample = pFileNameList[i * filesPerThread:]
                buffer_sample = pBuffer[i * filesPerThread:]

            queue[i] = Queue()
            process[i] = Process(target=writeSampleFiles, kwargs=dict(pFileNameList=file_list_sample,pBuffer=buffer_sample,pQueue=queue[i]))
            process[i].start()

        while not all_data_collected:
            for i in range(pThreads):
                if queue[i] is not None and not queue[i].empty():
                    _ = queue[i].get()
                    process[i].join()
                    process[i].terminate()
                    process[i] = None

            all_data_collected = True

            for i in range(pThreads):
                if process[i] is not None:
                    all_data_collected = False
            time.sleep(1)

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

with open(args.top_cell_barcode_list, 'r') as f:
    line = f.read().strip().split('\n')

barcode_set= set()
for i in range(len(line)):
    tmp = line[i].split()
    barcode_set.add(tmp[0])

# read and split fastq
count,line_count,cell_counter = 0,0,0
cell_index = {}

with gzip.open(args.R1,'rt') as f1,gzip.open(args.R2,'rt') as f2:
    lines_out_buffer,file_writer=[],[]
    pBufferSize, pThreads = int(args.buffer_size), int(args.thread)

    for line1,line2 in zip(f1,f2):
        if count % 4 == 0:
            # if line1.strip()[1:29] != line2.strip()[1:29]:
            #     print('readID from R1 not equal to readID from R2')
            # read_id = line1.strip()[1:29]
            read_id = line1.strip().split()[0]
        if count % 4 == 1:
            seq = line1.strip()+line2.strip()
            sample_index1, i7index = seq[200:206], seq[210:218]
        if count % 4 == 3:
            qual = line1.strip()+line2.strip()
            R1_insert_qual, R2_insert_qual, R1_insert, R2_insert = qual[3:100], qual[103:200], seq[3:100], seq[103:200]
            read_id_w_tag = read_id + '_' + sample_index1 + '_' + i7index
            barcode_tag = sample_index1 + '_' + i7index
            if barcode_tag in barcode_set:
                if barcode_tag in cell_index:
                    lines_out_buffer[cell_index[barcode_tag]][0].extend([read_id_w_tag+ '\n', R1_insert+ '\n', '+\n', R1_insert_qual+ '\n'])
                    lines_out_buffer[cell_index[barcode_tag]][1].extend([read_id_w_tag+ '\n', R2_insert+ '\n', '+\n', R2_insert_qual+ '\n'])
                else:
                    cell_index[barcode_tag] = cell_counter
                    file_writer.append([os.path.abspath(args.out_dir) + '/' + barcode_tag + '_R1.fq.gz',
                                        os.path.abspath(args.out_dir) + '/' + barcode_tag + '_R2.fq.gz'])
                    lines_out_buffer.append([[read_id_w_tag+ '\n', R1_insert+ '\n', '+\n', R1_insert_qual+ '\n'],
                                             [read_id_w_tag+ '\n', R2_insert+ '\n', '+\n', R2_insert_qual+ '\n']])
                    cell_counter += 1

                demultiplex_count[meta_dict[index1_sample_dict[sample_index1],i7_indexNum_dict[i7index]]] += 1
        count += 1
        line_count += 4

        if line_count % pBufferSize == 0 or line1 is False:
            buffered_elements = 0
            for lines_buffered in lines_out_buffer:
                buffered_elements += len(lines_buffered[0])
                buffered_elements += len(lines_buffered[1])
            if buffered_elements > pBufferSize or line1 is False:
                handleCompressingMulticore(file_writer, lines_out_buffer, pThreads)
                lines_out_buffer, file_writer=[],[]
                line_count = 0
                cell_index = {}
                cell_counter = 0

    if lines_out_buffer is not None and len(lines_out_buffer) > 0:
        handleCompressingMulticore(file_writer, lines_out_buffer, pThreads)
        lines_out_buffer, file_writer = [], []
        cell_counter = 0
        cell_index = {}

total_demultiplex_count = 0
for k in demultiplex_count:
    print(k)
    print(demultiplex_count[k])
    total_demultiplex_count += demultiplex_count[k]

print("total demultiplex count")
print(total_demultiplex_count)
print("total read count")
print(int(count/4))
