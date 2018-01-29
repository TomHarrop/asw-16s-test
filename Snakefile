#!/usr/bin/env python3



#############
# FUNCTIONS #
#############

###########
# GLOBALS #
###########

# primers
head = "GCTATGCGCGAGCTGC"
primer_f = "CMGGATTAGATACCCKGG"     # 799F-mod3 18nt
primer_r = "ACGGGCGGTGTGTRC"        # 1392R

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
        dynamic('output/matched_merged/{individual}_r1.fq'),
        dynamic('output/matched_merged/{individual}_r2.fq')

# merged the antisense and sense reads
rule merge_sense_antisense:
    input:
        r1_r1 = 'output/matched_sense/{individual}_r1.fq',
        r2_r1 = 'output/matched_antisense_renamed/{individual}_r1.fq',
        r1_r2 = 'output/matched_sense/{individual}_r2.fq',
        r2_r2 = 'output/matched_antisense_renamed/{individual}_r1.fq'
    output:
        r1 = 'output/matched_merged/{individual}_r1.fq',
        r2 = 'output/matched_merged/{individual}_r2.fq'
    threads:
        1
    script:
        'src/get_unique_reads.py'

# reverse complement the antisense matches
rule rename_antisense_reads:
    input:
        r1 = 'output/matched_antisense/{individual}_r1.fq',
        r2 = 'output/matched_antisense/{individual}_r2.fq',
    output:
        r1 = temp('output/matched_antisense_renamed/{individual}_r1.fq'),
        r2 = temp('output/matched_antisense_renamed/{individual}_r2.fq')
    threads:
        1
    script:
        'src/rename_antisense_reads.py'

# match the primers
rule match_sense_primers:
    input:
        r1 = 'output/matched_head/{individual}_r1.fq',
        r2 = 'output/matched_head/{individual}_r2.fq'
    output:
        r1 = temp('output/matched_sense/{individual}_r1.fq'),
        r2 = temp('output/matched_sense/{individual}_r2.fq')
    params:
        r1d = 'output/matched_sense/{individual}_r1.fq.discards',
        r2d = 'output/matched_sense/{individual}_r2.fq.discards'
    threads:
        1
    log:
        'output/logs/cutadapt-match-sense_{individual}.log'
    shell:
        'cutadapt '
        '-g 799F={primer_f} '
        '-G 1392R={primer_r} '
        '--minimum-length=207 '
        '--error-rate=0 '
        '--no-indels '
        '--pair-filter=any '
        '-o {output.r1} '
        '-p {output.r2} '
        '--untrimmed-output={params.r1d} '
        '--untrimmed-paired-output={params.r2d} '
        '{input.r1} {input.r2} '
        '&> {log}'

rule match_antisense_primers:
    input:
        r1 = 'output/matched_head/{individual}_r1.fq',
        r2 = 'output/matched_head/{individual}_r2.fq'
    output:
        r1 = temp('output/matched_antisense/{individual}_r1.fq'),
        r2 = temp('output/matched_antisense/{individual}_r2.fq')
    params:
        r1d = 'output/matched_antisense/{individual}_r1.fq.discards',
        r2d = 'output/matched_antisense/{individual}_r2.fq.discards'
    threads:
        1
    log:
        'output/logs/cutadapt-match-antisense_{individual}.log'
    shell:
        'cutadapt '
        '-g 1392R={primer_r} '
        '-G 799F={primer_f} '        
        '--minimum-length=207 '
        '--error-rate=0 '
        '--no-indels '
        '--pair-filter=any '
        '-o {output.r1} '
        '-p {output.r2} '
        '--untrimmed-output={params.r1d} '
        '--untrimmed-paired-output={params.r2d} '
        '{input.r1} {input.r2} '
        '&> {log}'

# match the HEAD sequence
rule match_head:
    input:
        r1 = 'output/demuxed_merged/{individual}_r1.fq',
        r2 = 'output/demuxed_merged/{individual}_r2.fq'
    output:
        r1 = temp('output/matched_head/{individual}_r1.fq'),
        r2 = temp('output/matched_head/{individual}_r2.fq')
    params:
        r1d = 'output/matched_head/{individual}_r1.fq.discards',
        r2d = 'output/matched_head/{individual}_r2.fq.discards'
    threads:
        1
    log:
        'output/logs/cutadapt-match-head_{individual}.log'
    shell:
        'cutadapt '
        '-g HEAD={head} '
        '-G HEAD={head} '
        '--minimum-length=225 '
        '--error-rate=0 '
        '--no-indels '
        '--pair-filter=any '
        '-o {output.r1} '
        '-p {output.r2} '
        '--untrimmed-output={params.r1d} '
        '--untrimmed-paired-output={params.r2d} '
        '{input.r1} {input.r2} '
        '&> {log}'

# merge and sort cudatapt output
rule merge_demuxed_reads:
    input:
        r1_r1 = 'output/cutadapt_demux_r1/{individual}_r1.fq',
        r2_r1 = 'output/cutadapt_demux_r2/{individual}_r1.fq',
        r1_r2 = 'output/cutadapt_demux_r1/{individual}_r2.fq',
        r2_r2 = 'output/cutadapt_demux_r2/{individual}_r2.fq'
    output:
        r1 = temp('output/demuxed_merged/{individual}_r1.fq'),
        r2 = temp('output/demuxed_merged/{individual}_r2.fq')
    threads:
        1
    script:
        'src/get_unique_reads.py'


# demux with cutadapt
rule demux_r1:
    input:
        key = 'output/config_files/cutadapt_barcodes.fasta',
        r1 = 'output/raw_reads/r1.fq',
        r2 = 'output/raw_reads/r2.fq'
    output:
        dynamic('output/cutadapt_demux_r1/{individual}_r1.fq'),
        dynamic('output/cutadapt_demux_r1/{individual}_r2.fq')
    log:
        'output/logs/cutadapt_demux_r1.log'
    threads:
        1
    shell:
        'cutadapt '
        '-g file:{input.key} '
        '--error-rate=0 '
        '--no-indels '
        '--no-trim '
        '--max-n=0 '
        '--pair-filter=any '
        '-o output/cutadapt_demux_r1/{{name}}_r1.fq '
        '-p output/cutadapt_demux_r1/{{name}}_r2.fq '
        '--untrimmed-output='
        'output/cutadapt_demux_r1/untrimmed_r1.fq.discards '
        '--untrimmed-paired-output='
        'output/cutadapt_demux_r1/untrimmed_r2.fq.discards '
        '{input.r1} {input.r2} '
        '&> {log}'

rule demux_r2:
    input:
        key = 'output/config_files/cutadapt_barcodes.fasta',
        r1 = 'output/raw_reads/r1.fq',
        r2 = 'output/raw_reads/r2.fq'
    output:
        dynamic('output/cutadapt_demux_r2/{individual}_r1.fq'),
        dynamic('output/cutadapt_demux_r2/{individual}_r2.fq')
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
        '-o output/cutadapt_demux_r2/{{name}}_r2.fq '
        '-p output/cutadapt_demux_r2/{{name}}_r1.fq '
        '--untrimmed-output='
        'output/cutadapt_demux_r2/untrimmed_r2.fq.discards '
        '--untrimmed-paired-output='
        'output/cutadapt_demux_r2/untrimmed_r1.fq.discards '
        '{input.r2} {input.r1} '
        '&> {log}'


# expand the reads to avoid pigz
rule zcat:
    input:
        'data/DDP02116-W/TH1_{r}.fq.gz'
    output:
        temp('output/raw_reads/r{r}.fq')
    threads:
        1
    shell:
        'zcat {input} > {output}'


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