#!/usr/bin/env python3

from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord

# This is quick and dirty. Split at 210 b because that's what I did during
# demuxing. Would be cleaner to take the ID and get the original reads


#############
# FUNCTIONS #
#############

# read splitter
def split_read(record, split_at):
    return record
    # my_spacer = Seq.Seq("N" * 205)
    # my_prefix = record[:split_at]
    # my_suffix = record[split_at:]
    # return my_prefix + my_spacer + my_suffix

###########
# GLOBALS #
###########

in_fasta = snakemake.input['fasta']
out_fasta = snakemake.output['fasta']

########
# MAIN #
########

records = [x for x in SeqIO.parse(in_fasta, 'fasta')]
split_records = [split_read(x, 210) for x in records]
SeqIO.write(split_records, out_fasta, 'fasta')
