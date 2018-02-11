#!/usr/bin/env python3

import csv
from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord


#############
# FUNCTIONS #
#############

def create_contig(read_id, fwd_reads, rev_reads):
    '''
    Use read_id to get a fwd_read from fwd_reads and a rev_read from rev_reads.
    Concatenate the fwd_read with the reverse complement of rev_read and return
    with and without a spacer.
    '''
    my_spacer = Seq.Seq('N')
    my_fwd_read = fwd_reads[read_id].seq
    my_rev_read = rev_reads[read_id].reverse_complement().seq
    return {
        'with_spacer': my_fwd_read + my_spacer + my_rev_read,
        'no_spacer':   my_fwd_read + my_rev_read}


def seqrec_wrapper(dna_sequence, record_id):
    return SeqRecord.SeqRecord(
        dna_sequence,
        id=record_id,
        name='',
        description='')

###########
# GLOBALS #
###########

r1 = snakemake.input['r1']
r2 = snakemake.input['r2']
individual = snakemake.params['individual']
fa_out = snakemake.output['fa']
spaced_fa_out = snakemake.output['spaced_fa']
out_key = snakemake.output['key']

########
# MAIN #
########

# get all identifiers and enumerate
r1_reads = SeqIO.to_dict(SeqIO.parse(r1, 'fastq-sanger'))
r2_reads = SeqIO.to_dict(SeqIO.parse(r2, 'fastq-sanger'))
id_to_key = {y: x for x, y in enumerate(r1_reads.keys())}

# generate records
spaced_records = []
unspaced_records = []
output_ids = {}
for rec_id in id_to_key.keys():
    my_contigs = create_contig(rec_id, r1_reads, r2_reads)
    my_rec_no = id_to_key[rec_id] + 1
    my_id = '{0}|{1}'.format(individual, str(my_rec_no))
    spaced_records.append(
        seqrec_wrapper(my_contigs['with_spacer'], my_id))
    unspaced_records.append(
        seqrec_wrapper(my_contigs['no_spacer'], my_id))
    output_ids[rec_id] = my_id

# write output
SeqIO.write(spaced_records, spaced_fa_out, 'fasta')
SeqIO.write(unspaced_records, fa_out, 'fasta')
with open(out_key, 'wt') as f:
    my_writer = csv.writer(f)
    my_writer.writerow(['read_id', 'contig_name'])
    my_writer.writerows(output_ids.items())
