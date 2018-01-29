#!/usr/bin/env python3

from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord



#############
# FUNCTIONS #
#############


# secrecord constructor
def create_contig(read_id, fwd_reads, rev_reads):
    my_spacer = Seq.Seq("N" * 205)
    my_fwd_read = fwd_reads[read_id].seq
    my_rev_read = rev_reads[read_id].reverse_complement().seq
    return my_fwd_read + my_spacer + my_rev_read


###########
# GLOBALS #
###########

r1 = snakemake.input['r1']
r2 = snakemake.input['r2']
individual = snakemake.params['individual']
out_file = snakemake.output['fa']

########
# MAIN #
########

# get all identifiers and enumerate
r1_reads = SeqIO.to_dict(SeqIO.parse(r1, 'fastq-sanger'))
r2_reads = SeqIO.to_dict(SeqIO.parse(r2, 'fastq-sanger'))
id_to_key = {y: x for x, y in enumerate(r1_reads.keys())}

# generate records
output_records = []
for rec_id in id_to_key.keys():
    my_seq = create_contig(rec_id, r1_reads, r2_reads)
    my_rec_no = id_to_key[rec_id] + 1
    my_id = '{0}|{1}'.format(individual, str(my_rec_no))
    output_records.append(SeqRecord.SeqRecord(
        my_seq,
        id=my_id,
        name='',
        description=''))

# write output
SeqIO.write(output_records, out_file, 'fasta')

