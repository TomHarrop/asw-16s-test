#!/usr/bin/env Rscript

library(data.table)
library(DESeq2)
library(ggplot2)
library(phyloseq)

# files
count_file <- "output/gutfilter/count_table.txt"
taxonomy_file <- "output/annotate_otus/keptotus.seed_v128.wang.taxonomy"

# generate OTU table
count_table <- fread(count_file)
count_matrix <- as.matrix(data.frame(count_table, row.names = "otu_id"))
otu <- otu_table(count_matrix, taxa_are_rows = TRUE)

# generate taxonomy matrix
SplitTaxonomy <- function(x) {
    no_nums <- gsub("\\([[:digit:]]+\\)", "", x)
    strsplit(no_nums, ";")
}
tax_levels <- c("Domain",
                "Phylum",
                "Class",
                "Order",
                "Family",
                "Genus")

taxonomy_table <- fread(taxonomy_file,
                        header = FALSE,
                        col.names = c("otu_id", "taxonomy"))
taxonomy_by_otu <- taxonomy_table[, transpose(SplitTaxonomy(taxonomy)),
                                  by = otu_id]
setnames(taxonomy_by_otu, paste0("V", c(1:6)), tax_levels)
tax_matrix <- as.matrix(data.frame(taxonomy_by_otu, row.names = "otu_id"))
tax <- tax_table(tax_matrix)

# generate sample data
sample_table <- data.table(samplename = colnames(count_matrix))
sample_table[, c("population", "individual") := tstrsplit(samplename, "_")]
sd <- sample_data(data.frame(sample_table, row.names = "samplename"))
sd$samplename <- rownames(sd)

# generate the phyloseq object
physeq <- phyloseq(otu, tax, sd)

# convert to DESeq2
dds <- phyloseq_to_deseq2(physeq, ~ population)

# run a wald test
dds_wald <- DESeq(dds)
res_wald <- results(dds_wald,
                    contrast = c("population", "Lincoln", "Invermay"),
                    lfcThreshold = log(1.5, 2))
subset(res_wald, padj < 0.1)
res_wald[order(res_wald$padj), ]
tax['Lincoln_17|3638']
plotCounts(dds_wald, "Lincoln_17|3638", intgroup = "population")

# run a likelihood ratio test on the populations
dds <- DESeq(dds, test = "LRT", reduced = ~ 1, fitType = "local")
res <- results(dds, alpha = 0.1, cooksCutoff = FALSE)
results_table <- data.table(data.frame(res), keep.rownames = TRUE)
setnames(results_table, "rn", "otu_id")


DESeq2::plotDispEsts(dds)
DESeq2::plotMA(dds)
vst_blind <- varianceStabilizingTransformation(
    dds, 
    blind = TRUE, fitType = "mean")


pca <- prcomp(t(assay(vst_blind)))

pca_table <- data.table(pca$x, keep.rownames = TRUE)
pca_table[, c("population", "individual") := tstrsplit(rn, "_")]
pca_pd <- melt(pca_table, id.vars = c("rn", "population", "individual"))
pca_pd[, dimension := as.numeric(gsub("[^[:digit:]]", "", variable))]

ggplot(pca_pd[dimension < 10],
       aes(x = population, y = value)) +
    facet_wrap(~ variable) +
    geom_point()

estimate_richness(physeq, measures = "Chao1")
sn_hc <- hclust(distance(physeq, method = "wunifrac"),
                method = "average")
plot(sn_hc)


dist_mat <- distance(physeq, method = "wunifrac")

