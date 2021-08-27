import json

import constraints

from multiprocessing import Pool
from itertools import zip_longest

import dinopy

fqr_1 = dinopy.FastqReader(str(snakemake.input[0]))
if(snakemake.params.paired):
    fqr_2 = dinopy.FastqReader(str(snakemake.input[1]))
else:
    fqr_2 = list()

filtered_1 = dinopy.FastqWriter(str(snakemake.output[0]))
if not snakemake.params.paired:
    filtered_out_1 = dinopy.FastqWriter(str(snakemake.output[1]))
else:
    filtered_2 = dinopy.FastqWriter(str(snakemake.output[1]))
    filtered_out_1 = dinopy.FastqWriter(str(snakemake.output[2]))
    filtered_out_2 = dinopy.FastqWriter(str(snakemake.output[3]))


filtered_1.open()
filtered_out_1.open()
if snakemake.params.paired:
    filtered_2.open()
    filtered_out_2.open()


def filter(reads):
    sequences = []
    sequences.append(reads[0][0].decode())
    if snakemake.params.paired:
        sequences.append(reads[1][0].decode())
    for i,sequence in enumerate(sequences):
        if snakemake.params.primer_length > 0:
            sequence = sequence[snakemake.params.primer_length:snakemake.params.sequence_length+snakemake.params.primer_length]
        for constraint in constraints.constraints:
            if (constraint.__name__ in snakemake.params.constraints['constraints']):
                if (constraint.__name__ == 'undesired_subsequences'):
                    failed = constraint(sequence, snakemake.input.subsequences_file)
                else:
                    failed = constraint(sequence, **snakemake.params.constraints[constraint.__name__])
                if failed:
                    return reads, failed, constraint.__name__, i+1
    return reads, False, '', None

def reads(reads):
    if isinstance(reads, list):
        return
    else:
        for read in reads.reads(quality_values=True):
            yield [read.sequence, read.name, read.quality]

pool = Pool(snakemake.threads)
results = pool.imap_unordered(filter, zip_longest(reads(fqr_1), reads(fqr_2)), 10000)
countdict = dict()
countdict[1] = { i : 0 for i in snakemake.params.constraints['constraints']}
countdict[1]['not_filtered'] = 0
countdict[1]['total'] = 0

if snakemake.params.paired:
    countdict[2] = {i: 0 for i in snakemake.params.constraints['constraints']}
    countdict[2]['not_filtered'] = 0
    countdict[2]['total'] = 0

for reads, failed, origin, r in results:
    countdict[1]['total'] += 1
    if snakemake.params.paired:
        countdict[2]['total'] += 1
    if failed:
        countdict[r][origin] += 1
        filtered_out_1.write(reads[0][0], reads[0][1], reads[0][2])
        if snakemake.params.paired:
            filtered_out_2.write(reads[1][0], reads[1][1], reads[1][2])
    else:
        countdict[1]['not_filtered'] += 1
        filtered_1.write(reads[0][0], reads[0][1], reads[0][2])
        if snakemake.params.paired:
            countdict[2]['not_filtered'] += 1
            filtered_2.write(reads[1][0], reads[1][1], reads[1][2])

with open(str(snakemake.log), 'w') as file:
    json.dump(countdict, file)
pool.close()  # 'TERM'
pool.join()   # 'KILL'

filtered_1.close()
filtered_out_1.close()
if snakemake.params.paired:
    filtered_2.close()
    filtered_out_2.close()
