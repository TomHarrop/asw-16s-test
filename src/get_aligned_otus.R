#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

###########
# GLOBALS #
###########

align_report_file <- snakemake@input[["report"]]
aligned_otu_file <- snakemake@output[["aligned_otus"]]
plot_file <- snakemake@output[["plot_file"]]
log_file <- snakemake@log[["log"]]

# dev
# align_report_file <- "output/091_annotate_otus/keptotus.align.report"
# aligned_otu_file <- "test_otus.txt"
# plot_file <- "test.pdf"

########
# MAIN #
########

# set log
log <- file(log_file, open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

# read table
align_report <- fread(align_report_file)

# calculate cutoffs
minstart <- align_report[, quantile(TemplateStart, 0.05)]
maxend <- align_report[, quantile(TemplateEnd, 0.95)]
xl <- c(
    align_report[, quantile(TemplateStart, 0.01)],
    align_report[, quantile(TemplateEnd, 0.99)]
)

# extract read names
aligned_reads <- align_report[TemplateStart >= minstart &
                                  TemplateEnd <= maxend,
                              .(unique(QueryName))]

# plot distributions
pd <- melt(align_report,
           id.vars = "QueryName",
           measure.vars = c("TemplateStart", "TemplateEnd"))

gt <- paste(dim(aligned_reads)[1],
      "OTUs aligned between template position",
      round(minstart, 0),
      "and",
      round(maxend, 0))

gp <- ggplot(pd, aes(x = value)) +
    ggtitle(gt) +
    xlab("Template position") + ylab("Count") +
    xlim(xl) +
    geom_vline(xintercept = c(minstart - 1, maxend + 1)) +
    facet_grid(variable ~ .) +
    geom_histogram(bins = 200)

# write output
fwrite(aligned_reads,
       aligned_otu_file)

ggsave(filename = plot_file,
       plot = gp,
       width = 10,
       height = 7.5,
       units = "in")

# write session info
sessionInfo()
