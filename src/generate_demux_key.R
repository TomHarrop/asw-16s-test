#!/usr/bin/env Rscript

library(data.table)
library(seqinr)

############
# FUNCTION #
############

ReverseCompString <- function(x) {
    toupper(
        paste0(
            rev(
                comp(
                    unlist(
                        strsplit(x,
                                 "")))),
            collapse = ""))
}

###########
# GLOBALS #
###########

barcodes_file <- snakemake@input[["barcodes_file"]]
key_file <- snakemake@input[["key_file"]]
key <- snakemake@output[["key"]]

#dev
# barcodes_file <- "data/bc_5nt_with_spacer.fasta"
# key_file <- "data/barcodes.csv"
# key <- "test.txt"

########
# MAIN #
########

barcodes <- fread(barcodes_file,
                  header = FALSE,
                  col.names = c("barcode_name", "sequence"))
sample_key <- fread(key_file)

# munge the barcodes
barcodes[, barcode_number := as.numeric(gsub("[^[:digit:]]+",
                                             "",
                                             barcode_name))]

# munge the sample key
sample_key[`Sample name` == "", `Sample name` := Sample]

# merge
sample_barcodes <- merge(
    barcodes[, .(bc = barcode_number,
                 sequence)],
    sample_key[, .(bc = Barcode,
                   sample = `Sample name`)])

# vector of names and barcodes
head_seq <- "GCTATGCGCGAGCTGC"
barcode_fasta <- sample_barcodes[, c(
    paste0(">", sample),
    paste0("^", sequence, head_seq),
    paste0(">", sample),
    paste0("^N", sequence, head_seq)),
    by = sample][, V1]

# write fasta
writeLines(barcode_fasta, key, sep = "\n")
