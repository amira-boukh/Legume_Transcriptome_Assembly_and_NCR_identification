#!/usr/bin/env Rscript

#Differential gene expression for Indigofera argentea nodules vs roots
#Requires Salmon quant files and a tx2gene mapping file
#Outputs DESeq2 and edgeR results tables and visualisations

library(tximport)
library(DESeq2)
library(tidyverse)
library(edgeR)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)
library(grid)

#Edit these variables before running

quant_table_file <- "DEGs/Quants.csv"       #columns: Sample, quant_file
tx2gene_file <- "DEGs/geneID.csv"           #two columns: transcript, gene
output_folder <- "DEGs/results_indigo"

#NCR/NCR-like genes to highlight (Optional)
ncr_genes <- c("Cluster-4043.56207","Cluster-3230.0","Cluster-4043.19630",
               "Cluster-4043.23593","Cluster-4043.48831","Cluster-4043.22256",
               "Cluster-4043.21323","Cluster-4043.23067","Cluster-4043.26814",
               "Cluster-4043.44242","Cluster-4043.30772","Cluster-4043.29463")
leghemoglobin_genes <- c("Cluster-4043.20145","Cluster-4043.23597","Cluster-4043.42100")
ahl_genes <- c("Cluster-4043.22527","Cluster-4043.22009")

#Load data

if (!file.exists(quant_table_file)) {
  stop("Cannot find quant table: ", quant_table_file)
}
if (!file.exists(tx2gene_file)) {
  stop("Cannot find tx2gene mapping: ", tx2gene_file)
}
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

indigo_samples <- read_delim(quant_table_file, show_col_types = FALSE)
tx2gene_map <- read_delim(tx2gene_file,
                          col_names = c("transcript", "gene"),
                          show_col_types = FALSE)

if (!all(c("Sample", "quant_file") %in% names(indigo_samples))) {
  stop("quant table must contain columns named 'Sample' and 'quant_file'")
}
if (!all(file.exists(indigo_samples$quant_file))) {
  missing <- indigo_samples$quant_file[!file.exists(indigo_samples$quant_file)]
  stop("Missing Salmon quant files:\n", paste(missing, collapse = "\n"))
}

txi <- tximport(files = indigo_samples$quant_file,
                type = "salmon",
                tx2gene = tx2gene_map,
                countsFromAbundance = "no")

colnames(txi$counts) <- indigo_samples$Sample

sampleTable <- data.frame(
  sampleName = colnames(txi$counts),
  condition = c("Root", "Root", "Root", "Nodules", "Nodules", "Nodules"),
  row.names = colnames(txi$counts)
)
sampleTable$condition <- factor(sampleTable$condition,
                                levels = c("Root", "Nodules"))

#DESeq2 workflow

dds <- DESeqDataSetFromTximport(txi = txi,
                                colData = sampleTable,
                                design = ~ condition)
dds <- dds[rowSums(counts(dds)) > 0, ]
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "Nodules", "Root"))
res_table <- as.data.frame(res) %>%
  rownames_to_column("gene_id") %>%
  arrange(padj, desc(abs(log2FoldChange)))

write_tsv(res_table, file.path(output_folder, "deseq2_all_genes.tsv"))

res_sig <- res_table %>%
  filter(!is.na(padj), padj <= 0.05, abs(log2FoldChange) > 1)
write_tsv(res_sig, file.path(output_folder, "deseq2_significant.tsv"))

#Optional shrinkage for visualisation (requires apeglm if available)
if (requireNamespace("apeglm", quietly = TRUE)) {
  resLFC <- lfcShrink(dds, coef = "condition_Nodules_vs_Root", type = "apeglm")
} else {
  resLFC <- res
}

png(file.path(output_folder, "MA_plot.png"), width = 1600, height = 1400, res = 300)
plotMA(res, ylim = c(-15, 15), main = "DESeq2 MA-plot")
abline(h = c(-1, 1), col = "red", lty = 2)
dev.off()

# edgeR quasi-likelihood (matches diCoExpress GLM step)

dge <- DGEList(counts = txi$counts)
keep_edge <- rowSums(cpm(dge) > 1) >= 3
dge <- dge[keep_edge, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
design_edge <- model.matrix(~ condition, data = sampleTable)
dge <- estimateDisp(dge, design_edge)
fit <- glmQLFit(dge, design_edge)
qlf <- glmQLFTest(fit, coef = 2)

edge_table <- topTags(qlf, n = Inf, sort.by = "none")$table %>%
  rownames_to_column("gene_id") %>%
  rename(log2FoldChange = logFC,
         pvalue = PValue,
         padj = FDR) %>%
  arrange(padj, desc(abs(log2FoldChange)))

write_tsv(edge_table, file.path(output_folder, "edgeR_all_genes.tsv"))

edge_sig <- edge_table %>%
  filter(!is.na(padj), padj <= 0.05, abs(log2FoldChange) > 1)
write_tsv(edge_sig, file.path(output_folder, "edgeR_significant.tsv"))

overlap_genes <- intersect(res_sig$gene_id, edge_sig$gene_id)
write_tsv(tibble(gene_id = overlap_genes),
          file.path(output_folder, "overlap_deseq2_edgeR.tsv"))

#Gene-level visualisations

plots_dir <- file.path(output_folder, "plots")
if (!dir.exists(plots_dir)) dir.create(plots_dir)

#Plot counts for NCR / marker genes
genes_to_plot <- unique(c(ncr_genes, leghemoglobin_genes, ahl_genes))
for (gene_id in genes_to_plot) {
  if (gene_id %in% rownames(dds)) {
    png(file.path(plots_dir, paste0("plotCounts_", gene_id, ".png")),
        width = 1400, height = 1200, res = 300)
    plotCounts(dds, gene = gene_id, intgroup = "condition",
               main = gene_id)
    dev.off()
  }
}

#Normalised expression heatmap
vsd <- vst(dds)
annot_col <- sampleTable %>%
  select(condition)

save_pheatmap_png <- function(x, filename, width = 6000, height = 4000, res = 600) {
  png(filename, width = width, height = height, res = res)
  grid.newpage()
  grid.draw(x$gtable)
  dev.off()
}

heatmap_genes <- intersect(ncr_genes, rownames(vsd))
if (length(heatmap_genes) > 1) {
  ph <- pheatmap(assay(vsd)[heatmap_genes, ],
                 cluster_rows = TRUE,
                 cluster_cols = FALSE,
                 show_rownames = TRUE,
                 legend = FALSE,
                 fontsize = 16,
                 annotation_col = annot_col,
                 annotation_colors = list(condition = c(Root = "#26b3ff",
                                                        Nodules = "#ffad73")),
                 color = brewer.pal(n = 9, name = "Greys"))
  save_pheatmap_png(ph, file.path(plots_dir, "heatmap_NCRs.png"))
}

#Volcano plots
res_clean <- res_table %>%
  filter(!is.na(padj)) %>%
  mutate(gene_type = case_when(
    log2FoldChange >= 1 & padj <= 0.05 ~ "up",
    log2FoldChange <= -1 & padj <= 0.05 ~ "down",
    TRUE ~ "ns"
  ))

cols <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey70")

volcano1 <- ggplot(res_clean,
                   aes(x = log2FoldChange, y = -log10(padj),
                       fill = gene_type, size = gene_type, alpha = gene_type)) +
  geom_point(shape = 21, colour = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  scale_fill_manual(values = cols) +
  scale_size_manual(values = c(up = 2, down = 2, ns = 1)) +
  scale_alpha_manual(values = c(up = 1, down = 1, ns = 0.5)) +
  labs(title = bquote("Differential expression in nodules vs roots of" ~ italic("I. argentea")),
       x = "log2 fold change (Nodules vs Root)",
       y = "-log10 adjusted p-value") +
  ylim(0, 500) +
  theme_bw() +
  theme(panel.grid = element_blank())

ggsave(file.path(plots_dir, "Volcano_overview.png"),
       plot = volcano1, width = 10, height = 7, dpi = 600)

#Annotated volcano with NCR / LegHB / AHL labels
res_clean$type <- "Other"
res_clean$type[res_clean$gene_id %in% ncr_genes] <- "NCR"
res_clean$type[res_clean$gene_id %in% leghemoglobin_genes] <- "LegHB"
res_clean$type[res_clean$gene_id %in% ahl_genes] <- "AHL"

highlight_genes <- res_clean %>%
  filter(type != "Other", padj <= 0.05)

volcano2 <- ggplot(res_clean,
                   aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = gene_type),
             alpha = 0.3, shape = 16, size = 1.2) +
  geom_point(data = highlight_genes,
             shape = 21, size = 2.5, fill = "firebrick", colour = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_label_repel(data = highlight_genes,
                   aes(label = type),
                   size = 3, force = 5, max.overlaps = Inf) +
  scale_colour_manual(values = cols) +
  labs(title = bquote("Highlighted NCR/AHL/LegHB genes in" ~ italic("I. argentea")),
       x = "log2 fold change (Nodules vs Root)",
       y = "-log10 adjusted p-value",
       colour = "Expression change") +
  ylim(0, 500) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.grid = element_blank(),
        legend.position = "right",
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16))

ggsave(file.path(plots_dir, "Volcano_highlighted.png"),
       plot = volcano2, width = 10, height = 7, dpi = 600)

#Combined summary

summary_lines <- c(
  paste("DESeq2 significant genes:", nrow(res_sig)),
  paste("edgeR significant genes:", nrow(edge_sig)),
  paste("Overlap:", length(overlap_genes))
)
writeLines(summary_lines, file.path(output_folder, "summary.txt"))


