## 1. demux by barcode

- don't trim at this stage, that will be done in step 2
    + `--no-trim`
- Either R1 or R2 **must** have the complete barcode **with spacer and HEAD** in sense orientation
    + remove unmatched reads: `--untrimmed-output`
    + for demultiplexing, we can only search for the bc in the forward read read: `--pair-filter=any`
- 0 mismatches, 0 indels, but add one N to front of adaptor sequence
    + `-e 0`
    + `--no-indels`
- demultiplex with cutadapt:
    + specify the barcode file: `-g file:barcodes.fasta`
    + `-o trimmed{name}.fq`
    + name will be replaced by the name from the barcode file
    + need to run once to search for the barcode in R1 and once to look in R2

## 1b. merge matched reads

- `src/get_unique_reads.py`

## 2. match HEAD sequence

- both reads **must** have head sequence
    + in this case `--pair-filter=any` will discard read pairs where either read doesn't match the adaptor to `--untrimmed-output`
- the head must be near the three prime
    + to account for the barcode + HEAD (up to 8 + 16 b), `--minimum-length=225`
- 0 mismatches, 0 indels
    + `-e 0`
    + `--no-indels`

## 3. demux by F and R primer
- Either (R1 contains F+ **AND** R2 contains R+) **OR** (R1 contains R+ **AND** R2 contains F+)
    + in this case `--pair-filter=any` will discard read pairs where either read doesn't match the primer to `--untrimmed-output`
    + run once for FR amplicons and once for RF ('antisense') amplicons
- the head must be near the three prime
    + to account for the barcode + HEAD + primer (up to 8 + 16 + 18 b), `--minimum-length=207`
- 0 mismatches, 0 indels
    + `-e 0`
    + `--no-indels`

## 3b. swap R1 and R2 in the matched antisense reads

- `src/rename_antisense_reads.py`

## 3c. merge matched reads

- `src/get_unique_reads.py`
