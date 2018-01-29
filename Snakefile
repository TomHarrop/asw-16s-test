#!/usr/bin/env python3



#############
# FUNCTIONS #
#############

###########
# GLOBALS #
###########

# primers
head = "GCTATGCGCGAGCTGC"
primer_f = "GCTATGCGCGAGCTGCCMGGATTAGATACCCKGG"     # HEAD-799F-mod3
primer_r = "GCTATGCGCGAGCTGCACGGGCGGTGTGTRC"        # HEAD-1392R

# files and folders
barcodes_file = 'data/bc_5nt_with_spacer.fasta'
key_file = 'data/barcodes.csv'
r1 = 'data/DDP02116-W/TH1_1.fq.gz'
r2 = 'data/DDP02116-W/TH1_2.fq.gz'

#########
# RULES #
#########

rule target:
    input:
        dynamic('output/demuxed_merged/{individual}_r1.fq'),
        dynamic('output/demuxed_merged/{individual}_r2.fq')


# merge and sort cudatapt output
rule merge_demuxed_reads:
    input:
        r1_r1 = 'output/cutadapt_demux_r1/{individual}_r1.fq.gz',
        r2_r1 = 'output/cutadapt_demux_r2/{individual}_r1.fq.gz',
        r1_r2 = 'output/cutadapt_demux_r1/{individual}_r2.fq.gz',
        r2_r2 = 'output/cutadapt_demux_r2/{individual}_r2.fq.gz'
    output:
        r1 = 'output/demuxed_merged/{individual}_r1.fq',
        r2 = 'output/demuxed_merged/{individual}_r2.fq'
    threads:
        1
    script:
        'src/get_unique_reads.py'


# demux with cutadapt
rule demux_r1:
    input:
        key = 'output/config_files/cutadapt_barcodes.fasta',
        r1 = r1,
        r2 = r2
    output:
        dynamic('output/cutadapt_demux_r1/{individual}_r1.fq.gz'),
        dynamic('output/cutadapt_demux_r1/{individual}_r2.fq.gz')
    log:
        'output/logs/cutadapt_demux_r1.log'
    threads:
        1
    shell:
        'cutadapt '
        '-g file:{input.key} '
        '--error-rate=0 '
        '--no-indels '
        '--no-trim '                # ONLY FOR TESTING
        '--max-n=0 '
        '--pair-filter=any '
        '-o output/cutadapt_demux_r1/{{name}}_r1.fq.gz '
        '-p output/cutadapt_demux_r1/{{name}}_r2.fq.gz '
        '--untrimmed-output='
        'output/cutadapt_demux_r1/untrimmed_r1.fq.gz.discards '
        '--untrimmed-paired-output='
        'output/cutadapt_demux_r1/untrimmed_r2.fq.gz.discards '
        '{input.r1} {input.r2} '
        '&> {log}'

rule demux_r2:
    input:
        key = 'output/config_files/cutadapt_barcodes.fasta',
        r1 = r1,
        r2 = r2
    output:
        dynamic('output/cutadapt_demux_r2/{individual}_r1.fq.gz'),
        dynamic('output/cutadapt_demux_r2/{individual}_r2.fq.gz')
    log:
        'output/logs/cutadapt_demux_r2.log'
    threads:
        1
    shell:
        'cutadapt '
        '-g file:{input.key} '
        '--error-rate=0 '
        '--no-indels '
        '--no-trim '                # ONLY FOR TESTING
        '--max-n=0 '
        '--pair-filter=any '
        '-o output/cutadapt_demux_r2/{{name}}_r2.fq.gz '
        '-p output/cutadapt_demux_r2/{{name}}_r1.fq.gz '
        '--untrimmed-output='
        'output/cutadapt_demux_r2/untrimmed_r2.fq.gz.discards '
        '--untrimmed-paired-output='
        'output/cutadapt_demux_r2/untrimmed_r1.fq.gz.discards '
        '{input.r2} {input.r1} '
        '&> {log}'


# generate a key file for demuxing
rule generate_demux_key:
    input:
        barcodes_file = barcodes_file,
        key_file = key_file
    threads:
        1
    output:
        key = 'output/config_files/cutadapt_barcodes.fasta'
    script:
        'src/generate_demux_key.R'