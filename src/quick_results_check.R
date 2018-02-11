#!/usr/bin/env Rscript

library(data.table)
library(DESeq2)
library(ggplot2)
library(phyloseq)
library(scales)

# files
count_file <- "output/082_gutfilter/count_table.txt"
taxonomy_file <- "output/091_annotate_otus/keptotus.seed_v128.wang.taxonomy"
tree_file <- "output/100_tree/keptotus.phylip.taxonomy.tre"

# generate OTU table
count_table <- fread(count_file)
sdc <- count_table[, names(count_table)[sapply(.SD, is.numeric)]]
count_table[, lapply(.SD, sum), .SDcols = sdc]
count_table[, c("ASW13_1", "ASW13_2", "ASW13_3") := NULL] # try without dummy sample
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

# generate tree
tree <- phyloseq::read_tree(tree_file)

# generate the phyloseq object
physeq <- phyloseq(otu, tax, sd, tree)
physeq <- subset_taxa(physeq, !grepl("unclassified", Family)) # try removing unclassified

# convert to DESeq2
dds <- phyloseq_to_deseq2(physeq, ~ population)

# run a wald test
dds_wald <- DESeq(dds)
res_wald <- results(dds_wald,
                    contrast = c("population", "Lincoln", "Invermay"),
                    lfcThreshold = log(1.5, 2))
subset(res_wald, padj < 0.1)
res_wald[order(res_wald$padj), ]

tax[rownames(subset(res_wald, padj < 0.05))]

tax['Invermay_15|24']
tax['Invermay_15|2']
plotCounts(dds_wald, "Invermay_15|24", intgroup = "population")

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

plotPCA(vst_blind, intgroup = "population")
pca <- prcomp(t(assay(vst_blind)))

pca_table <- data.table(pca$x, keep.rownames = TRUE)
pca_table[, c("population", "individual") := tstrsplit(rn, "_")]
pca_pd <- melt(pca_table, id.vars = c("rn", "population", "individual"))
pca_pd[, dimension := as.numeric(gsub("[^[:digit:]]", "", variable))]

ggplot(pca_pd[dimension < 10],
       aes(x = population, y = value, colour = individual)) +
    facet_wrap(~ variable) +
    geom_point(position = position_jitter(width = 0.1))

estimate_richness(physeq, measures = "Chao1")
sn_hc <- hclust(distance(physeq, method = "wunifrac"),
                method = "average")
plot(sn_hc)


dist_mat <- distance(physeq, method = "wunifrac")



counts_wide <- counts(dds, normalized = TRUE)
norm_counts <- melt(data.table(counts_wide, keep.rownames = TRUE),
                    id = 'rn',
                    variable.name = "Sample",
                    value.name = "normalized_counts")
setnames(norm_counts, "rn", "OTU")

norm_counts_tax <- merge(norm_counts,
                         taxonomy_by_otu,
                         by.x = "OTU",
                         by.y = "otu_id",
                         all.x = TRUE)

norm_counts_by_genus <- norm_counts_tax[, .(total_reads = sum(normalized_counts,
                                                              na.rm = TRUE)),
                                        by = .(Sample, Genus, Phylum)]
norm_counts_by_genus[, sum_by_genus := sum(total_reads), by = Genus]
norm_counts_by_genus[, sum_by_sample := sum(total_reads), by = Sample]
norm_counts_by_genus[, pct_reads := total_reads * 100 / sum_by_sample]
setorder(norm_counts_by_genus, Phylum, -sum_by_genus)
norm_counts_by_genus[, Genus := factor(Genus, levels = unique(Genus))]

plot_genus <- norm_counts_by_genus[!grepl("unclassified", Genus),
                                   any(total_reads > 10),
                                   by = Genus][
                                       V1 == TRUE, unique(Genus)]
plot_otu <- norm_counts_tax[Genus %in% as.character(plot_genus), unique(OTU)]

bp_pd <- norm_counts_by_genus[Genus %in% plot_genus]

ggplot(bp_pd[total_reads > 0],
             aes(x = Genus, y = total_reads, fill = Phylum)) +
    theme_minimal(base_size = 12) +
    #scale_fill_brewer(palette = "Set1") +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    
    theme(axis.text.x = element_text(angle = 45,
                                     vjust = 1,
                                     hjust = 1),
          strip.background = element_blank(),
          strip.text.x = element_blank()) +
    ylab("Abundance") +
    facet_grid(Sample ~ Phylum, scales = "free_x", space = "free_x") +
    geom_col(position = "dodge", width = 0.5, size = 0.5)


##################
# PHYLOSEQ PLOTS #
##################


# gapstats, no idea what is
sord <- ordinate(physeq, method = "PCoA", distance = "wunifrac")
gs <- gapstat_ord(sord, axes=1:4, K.max = 4)
phyloseq::plot_clusgap(gs)

# ordination etc.
phyloseq::plot_heatmap(physeq, method = "wunifrac")
phyloseq::plot_net(physeq, point_label = "samplename")
phyloseq::plot_ordination(physeq, sord, type="samples", color="population") +
    theme_minimal()

plot_bar(subset_taxa(physeq, !grepl("unclassified", Family)),
         fill = "Family",
         x = "samplename") +
    facet_wrap(~ population, scales = "free_x")


plot_heatmap(physeq, taxa.label = "Phylum")

plot_richness(physeq)
estimate_richness(physeq, measures = "Chao1")

plot_tree(physeq, color = "population",
          shape = "Phylum",
          size = "Abundance",
          label.tips = "Genus")

plot_tree(subset_taxa(physeq, !grepl("unclassified", Genus)),
          color = "population",
          size = "Abundance",
          label.tips = "Genus")

plot_tree(subset_taxa(physeq, !grepl("unclassified", Genus)),
          color = "Class",
          size = "Abundance",
          label.tips = "Genus")

plot_tree(subset_taxa(physeq, Order == "Enterobacteriales"),
          color = "population",
          shape = "Family",
          size = "Abundance",
          label.tips = "Genus")

plot_tree(subset_taxa(physeq, Class == "Gammaproteobacteria"),
          color = "population",
          shape = "Order",
          size = "Abundance",
          label.tips = "Genus")

plot_tree(subset_taxa(physeq, Class == "Alphaproteobacteria"),
          color = "population",
          shape = "Order",
          size = "Abundance",
          label.tips = "Genus")


plot_tree(subset_taxa(physeq, Order == "Xanthomonadales"),
          color = "individual",
          shape = "Family",
          size = "Abundance",
          label.tips = "Genus")


