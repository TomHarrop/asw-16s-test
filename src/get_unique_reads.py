#!/usr/bin/env python3

from Bio import SeqIO
import gzip

###########
# GLOBALS #
###########

# r1_file = 'output/cutadapt_demux_r1/ASW13_1.fq.gz'
# r2_file = 'output/cutadapt_demux_r2/ASW13_1.fq.gz'
# names_file = 'test.txt'

r1_file = snakemake.input['r1']
r2_file = snakemake.input['r2']
names_file = snakemake.output['kept_reads']

# read the record ids
with gzip.open(r1_file, 'rt') as handle:
    fastq = SeqIO.parse(handle, 'fastq-sanger')
    r1_names = set(x.name for x in fastq)

with gzip.open(r2_file, 'rt') as handle:
    fastq = SeqIO.parse(handle, 'fastq-sanger')
    r2_names = set(x.name for x in fastq)

# join the sets
all_names = '\n'.join(set(r1_names.union(r2_names)))

# write to file
with open(names_file, 'wt') as f:
    f.writelines(all_names)
