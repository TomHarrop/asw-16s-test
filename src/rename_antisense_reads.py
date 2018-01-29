#!/usr/bin/env python3

from Bio import SeqIO
import re


#############
# FUNCTIONS #
#############

def r1_to_r2(read_description):
    return re.sub("\s1", " 2", read_description)


def r2_to_r1(read_description):
    return re.sub("\s2", " 1", read_description)


###########
# GLOBALS #
###########

r1_in = snakemake.input['r1']
r2_in = snakemake.input['r2']
r1_out = snakemake.output['r1']
r2_out = snakemake.output['r2']

########
# MAIN #
########

# rename R1
r1_records = [x for x in SeqIO.parse(r1_in, 'fastq-sanger')]
for x in r1_records:
    x.description = r1_to_r2(x.description)

# rename R2
r2_records = [x for x in SeqIO.parse(r2_in, 'fastq-sanger')]
for x in r2_records:
    x.description = r2_to_r1(x.description)

# write output
SeqIO.write(
    r1_records,
    r2_out,
    'fastq-sanger')
SeqIO.write(
    r2_records,
    r1_out,
    'fastq-sanger')
