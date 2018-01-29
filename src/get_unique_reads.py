#!/usr/bin/env python3

from Bio import SeqIO


#############
# FUNCTIONS #
#############

def get_unique_reads(file_1, file_2):
    fastq = [x for x in SeqIO.parse(file_1, 'fastq-sanger')]
    fastq_names = [x.name for x in fastq]
    for x in SeqIO.parse(file_2, 'fastq-sanger'):
        if x.name not in fastq_names:
            fastq.append(x)
    return SeqIO.to_dict(fastq)

###########
# GLOBALS #
###########

r1_r1 = snakemake.input['r1_r1']
r2_r1 = snakemake.input['r2_r1']
r1_r2 = snakemake.input['r1_r2']
r2_r2 = snakemake.input['r2_r2']
r1_out = snakemake.output['r1']
r2_out = snakemake.output['r2']

########
# MAIN #
########

r1_unique = get_unique_reads(r1_r1, r2_r1)
r2_unique = get_unique_reads(r1_r2, r2_r2)

key_order = sorted(r1_unique.keys())

SeqIO.write((r1_unique[x] for x in key_order),
    r1_out,
    'fastq-sanger')
SeqIO.write((r2_unique[x] for x in key_order),
    r2_out,
    'fastq-sanger')
