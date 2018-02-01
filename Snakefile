#!/usr/bin/env python3

import pathlib

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
silva_align = pathlib.Path('data/silva/silva.seed_v128.align').resolve()
silva_tax = pathlib.Path('data/silva/silva.seed_v128.tax').resolve()


#########
# RULES #
#########

rule target:
    input:
        'output/091_annotate_otus/keptotus.seed_v128.wang.taxonomy'

# 09 put Ns in the OTUs and annotate with mothur
rule annotate_otus:
    input:
        fasta = 'output/091_annotate_otus/keptotus.fasta'
    output:
        tax = 'output/091_annotate_otus/keptotus.seed_v128.wang.taxonomy'
    params:
        wd = 'output/091_annotate_otus',
        align = silva_align,
        tax = silva_tax
    log:
        'output/logs/mothur_annotate-otus.log'
    threads:
        20
    shell:
        'bash -c \''
        'cd {params.wd} || exit 1 ; '
        'mothur "'
        '#classify.seqs(fasta=keptotus.fasta, '
        'template={params.align}, '
        'taxonomy={params.tax}, '
        'processors={threads})" '
        '\' &> {log}'

rule insert_ns:
    input:
        fasta = 'output/082_gutfilter/keptotus.fasta'
    output:
        fasta = 'output/091_annotate_otus/keptotus.fasta'
    threads:
        1
    script:
        'src/split_otus_for_mothur.py'  

# 082 run gutfilter and extract reads
rule gutfilter_reads:
    input:
        kept_otus = 'output/082_gutfilter/kept_otus.txt',
        fasta = 'output/081_swarm_reformatted/precluster.fasta'
    output:
        'output/082_gutfilter/keptotus.fasta'
    log:
        'output/logs/filterbyname.log'
    threads:
        1
    shell:
        'filterbyname.sh '
        'in={input.fasta} '
        'out={output} '
        'names={input.kept_otus} '
        'include=t '
        '2> {log}'

rule gutfilter:
    input:
        precluster_names = 'output/081_swarm_reformatted/long_table.tab'
    output:
        kept_otus = 'output/082_gutfilter/kept_otus.txt',
        count_table = 'output/082_gutfilter/count_table.txt',
        abundance_table = 'output/082_gutfilter/abundance_table.txt',
        filter_file = 'output/082_gutfilter/filter_table.txt'
    threads:
        1
    script:
        'src/gut_filter.R'

# 081 reformat for gutfilter / R
rule reformat_swarm_output:
    input:
        mothur_names = 'output/06_dereplicate/all.names',
        swarm_results = 'output/07_swarm/all.unique.swarm',
        dereplicated_fasta = 'output/06_dereplicate/all.unique.fasta'
    output:
        names = 'output/081_swarm_reformatted/precluster.names',
        fasta = 'output/081_swarm_reformatted/precluster.fasta',
        long_table = 'output/081_swarm_reformatted/long_table.tab'
    threads:
        1
    script:
        'src/reformat_swarm_output.py'

# 072 Cluster using swarm:
rule swarm:
    input:
        'output/07_swarm/all.unique.fasta'
    output:
        'output/07_swarm/all.unique.swarm'
    threads:
        20
    log:
        'output/logs/swarm.log'
    shell:
        'swarm '
        '-t {threads} '
        '-d 6 '
        # '--fastidious '
        # '-d 1 '
        '-l {log} '
        '-o {output} '
        '{input}'

# 071 Reformat mothur’s output for swarm:
rule reformat_for_swarm:
    input:
        fa = 'output/06_dereplicate/all.unique.fasta',
    output:
        fa = 'output/07_swarm/all.unique.fasta',
    params:
        names = 'output/06_dereplicate/all.names'
    threads:
        1
    script:
        'src/reheader_for_swarm.py'

# 06 Dereplicate using mothur
rule dereplicate:
    input:
        'output/054_joined_reads/all.fasta'
    output:
        fa = temp('output/06_dereplicate/all.fasta'),
        derep = 'output/06_dereplicate/all.unique.fasta',
        names = 'output/06_dereplicate/all.names'
    threads:
        1
    params:
        wd = 'output/06_dereplicate'
    log:
        'output/logs/dereplicate.log'
    shell:
        'cp {input} {output.fa} ; '
        'bash -c \''
        'cd {params.wd} || exit 1 ; '
        'mothur "'
        '#unique.seqs(fasta=all.fasta)" '
        '\' &> {log}'

# 054 join read files
rule join_reads:
    input:
        dynamic('output/053_renamed_contigs/{individual}.fa')
    output:
        'output/054_joined_reads/all.fasta'
    shell:
        'cat {input} > {output}'

# 053 make contigs and rename reads
rule contig_and_rename:
    input:
        r1 = 'output/052_truncated/{individual}_r1.fq',
        r2 = 'output/052_truncated/{individual}_r2.fq'        
    output:
        fa = 'output/053_renamed_contigs/{individual}.fa',
        key = 'output/053_renamed_contigs/{individual}_readkey.txt'
    params:
        individual = '{individual}'
    threads:
        1 
    script:
        'src/contig_and_rename.py'

# 05b truncate the reads to 200b
rule truncate:
    input:
        r1 = 'output/051_matched_merged/{individual}_r1.fq',
        r2 = 'output/051_matched_merged/{individual}_r2.fq'
    output:
        r1 = 'output/052_truncated/{individual}_r1.fq',
        r2 = 'output/052_truncated/{individual}_r2.fq'
    threads:
        20
    log:
        'output/logs/truncate_{individual}.log'
    shell:
        'bbduk.sh '
        'threads={threads} '
        'in={input.r1} '
        'in2={input.r2} '
        'qtrim=r '
        'trimq=10 '
        'maxns=0 '
        'literal='
        'AAAAAAAAAA,'
        'CCCCCCCCCC,'
        'GGGGGGGGGG,'
        'TTTTTTTTTT '
        'maskmiddle=f '
        'out={output.r1} '
        'out2={output.r2} '
        'forcetrimright=179 '
        'minlength=180 '
        '2> {log}'

# 05a merge the antisense and sense reads
rule merge_sense_antisense:
    input:
        r1_r1 = 'output/041_matched_sense/{individual}_r1.fq',
        r2_r1 = 'output/042_matched_antisense_renamed/{individual}_r1.fq',
        r1_r2 = 'output/041_matched_sense/{individual}_r2.fq',
        r2_r2 = 'output/042_matched_antisense_renamed/{individual}_r1.fq'
    output:
        r1 = 'output/051_matched_merged/{individual}_r1.fq',
        r2 = 'output/051_matched_merged/{individual}_r2.fq'
    threads:
        1
    script:
        'src/get_unique_reads.py'

# 04b reverse complement the antisense matches
rule rename_antisense_reads:
    input:
        r1 = 'output/041_matched_antisense/{individual}_r1.fq',
        r2 = 'output/041_matched_antisense/{individual}_r2.fq',
    output:
        r1 = 'output/042_matched_antisense_renamed/{individual}_r1.fq',
        r2 = 'output/042_matched_antisense_renamed/{individual}_r2.fq'
    threads:
        1
    script:
        'src/rename_antisense_reads.py'

# 04 match the primers
rule match_sense_primers:
    input:
        r1 = 'output/03_matched_head/{individual}_r1.fq',
        r2 = 'output/03_matched_head/{individual}_r2.fq'
    output:
        r1 = 'output/041_matched_sense/{individual}_r1.fq',
        r2 = 'output/041_matched_sense/{individual}_r2.fq'
    params:
        r1d = 'output/041_matched_sense/{individual}_r1.fq.discards',
        r2d = 'output/041_matched_sense/{individual}_r2.fq.discards'
    threads:
        1
    log:
        'output/logs/cutadapt-match-sense_{individual}.log'
    shell:
        'cutadapt '
        '-g 799F={primer_f} '
        '-G 1392R={primer_r} '
        '--minimum-length=207 '
        # '--error-rate=0.12 '
        # '--no-indels '
        '--pair-filter=any '
        '-o {output.r1} '
        '-p {output.r2} '
        '--untrimmed-output={params.r1d} '
        '--untrimmed-paired-output={params.r2d} '
        '{input.r1} {input.r2} '
        '&> {log}'

rule match_antisense_primers:
    input:
        r1 = 'output/03_matched_head/{individual}_r1.fq',
        r2 = 'output/03_matched_head/{individual}_r2.fq'
    output:
        r1 = 'output/041_matched_antisense/{individual}_r1.fq',
        r2 = 'output/041_matched_antisense/{individual}_r2.fq'
    params:
        r1d = 'output/041_matched_antisense/{individual}_r1.fq.discards',
        r2d = 'output/041_matched_antisense/{individual}_r2.fq.discards'
    threads:
        1
    log:
        'output/logs/cutadapt-match-antisense_{individual}.log'
    shell:
        'cutadapt '
        '-g 1392R={primer_r} '
        '-G 799F={primer_f} '        
        '--minimum-length=207 '
        # '--error-rate=0 '
        # '--no-indels '
        '--pair-filter=any '
        '-o {output.r1} '
        '-p {output.r2} '
        '--untrimmed-output={params.r1d} '
        '--untrimmed-paired-output={params.r2d} '
        '{input.r1} {input.r2} '
        '&> {log}'

# 03 match the HEAD sequence
rule match_head:
    input:
        r1 = 'output/022_demuxed_merged/{individual}_r1.fq',
        r2 = 'output/022_demuxed_merged/{individual}_r2.fq'
    output:
        r1 = 'output/03_matched_head/{individual}_r1.fq',
        r2 = 'output/03_matched_head/{individual}_r2.fq'
    params:
        r1d = 'output/03_matched_head/{individual}_r1.fq.discards',
        r2d = 'output/03_matched_head/{individual}_r2.fq.discards'
    threads:
        1
    log:
        'output/logs/cutadapt-match-head_{individual}.log'
    shell:
        'cutadapt '
        '-g HEAD={head} '
        '-G HEAD={head} '
        '--minimum-length=225 '
        # '--error-rate=0 '
        # '--no-indels '
        '--pair-filter=any '
        '-o {output.r1} '
        '-p {output.r2} '
        '--untrimmed-output={params.r1d} '
        '--untrimmed-paired-output={params.r2d} '
        '{input.r1} {input.r2} '
        '&> {log}'

# 02b merge and sort cudatapt output
rule merge_demuxed_reads:
    input:
        r1_r1 = 'output/021_cutadapt_demux_r1/{individual}_r1.fq',
        r2_r1 = 'output/021_cutadapt_demux_r2/{individual}_r1.fq',
        r1_r2 = 'output/021_cutadapt_demux_r1/{individual}_r2.fq',
        r2_r2 = 'output/021_cutadapt_demux_r2/{individual}_r2.fq'
    output:
        r1 = 'output/022_demuxed_merged/{individual}_r1.fq',
        r2 = 'output/022_demuxed_merged/{individual}_r2.fq'
    threads:
        1
    script:
        'src/get_unique_reads.py'


# 02 demux with cutadapt
rule demux_r1:
    input:
        key = 'output/01_config_files/cutadapt_barcodes.fasta',
        r1 = 'output/01_raw_reads/r1.fq',
        r2 = 'output/01_raw_reads/r2.fq'
    output:
        dynamic('output/021_cutadapt_demux_r1/{individual}_r1.fq'),
        dynamic('output/021_cutadapt_demux_r1/{individual}_r2.fq')
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
        '-o output/021_cutadapt_demux_r1/{{name}}_r1.fq '
        '-p output/021_cutadapt_demux_r1/{{name}}_r2.fq '
        '--untrimmed-output='
        'output/021_cutadapt_demux_r1/untrimmed_r1.fq.discards '
        '--untrimmed-paired-output='
        'output/021_cutadapt_demux_r1/untrimmed_r2.fq.discards '
        '{input.r1} {input.r2} '
        '&> {log}'

rule demux_r2:
    input:
        key = 'output/01_config_files/cutadapt_barcodes.fasta',
        r1 = 'output/01_raw_reads/r1.fq',
        r2 = 'output/01_raw_reads/r2.fq'
    output:
        dynamic('output/021_cutadapt_demux_r2/{individual}_r1.fq'),
        dynamic('output/021_cutadapt_demux_r2/{individual}_r2.fq')
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
        '-o output/021_cutadapt_demux_r2/{{name}}_r2.fq '
        '-p output/021_cutadapt_demux_r2/{{name}}_r1.fq '
        '--untrimmed-output='
        'output/021_cutadapt_demux_r2/untrimmed_r2.fq.discards '
        '--untrimmed-paired-output='
        'output/021_cutadapt_demux_r2/untrimmed_r1.fq.discards '
        '{input.r2} {input.r1} '
        '&> {log}'


# 01b expand the reads to avoid pigz
rule zcat:
    input:
        'data/DDP02116-W/TH1_{r}.fq.gz'
    output:
        temp('output/01_raw_reads/r{r}.fq')
    threads:
        1
    shell:
        'zcat {input} > {output}'


# 01a generate a key file for demuxing
rule generate_demux_key:
    input:
        barcodes_file = barcodes_file,
        key_file = key_file
    threads:
        1
    output:
        key = 'output/01_config_files/cutadapt_barcodes.fasta'
    script:
        'src/generate_demux_key.R'