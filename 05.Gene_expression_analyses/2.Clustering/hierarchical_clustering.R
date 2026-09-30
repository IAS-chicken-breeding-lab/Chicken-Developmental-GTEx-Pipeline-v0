#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggtree)
})

args <- commandArgs(trailingOnly = TRUE)

expr_file  <- args[1]
meta_file  <- args[2]
color_file <- args[3]
out_file   <- args[4]

# Read data
expr <- fread(expr_file, data.table = FALSE)
meta <- fread(meta_file, data.table = FALSE)
cols <- fread(color_file, data.table = FALSE)

colnames(expr)[1] <- "Gene_ID"

if ("Peroid" %in% colnames(meta))
  colnames(meta)[colnames(meta) == "Peroid"] <- "Period"

# Match samples
common <- intersect(colnames(expr)[-1], meta$BioSample)

if (length(common) == 0)
  stop("No matched samples.")

expr <- expr[, c("Gene_ID", common), drop = FALSE]
meta <- meta[match(common, meta$BioSample), , drop = FALSE]

# Expression matrix
mat <- as.matrix(expr[, -1])
storage.mode(mat) <- "numeric"
rownames(mat) <- expr$Gene_ID

# Log transformation
mat <- log2(mat + 0.25)

# Top 5000 highly variable genes
gene_sd <- apply(mat, 1, sd, na.rm = TRUE)
gene_sd[is.na(gene_sd)] <- 0

top_genes <- order(gene_sd, decreasing = TRUE)
top_genes <- top_genes[seq_len(min(5000, length(top_genes)))]

mat <- mat[top_genes, , drop = FALSE]

# Pearson correlation + hierarchical clustering
cor_mat <- cor(mat, method = "pearson", use = "pairwise.complete.obs")

tree <- hclust(
  as.dist(1 - cor_mat),
  method = "complete"
)

# Colors
tissue_col <- setNames(
  cols$color[cols$Category == "Tissue"],
  cols$sample[cols$Category == "Tissue"]
)

period_col <- setNames(
  cols$color[cols$Category == "Period"],
  cols$sample[cols$Category == "Period"]
)

sex_col <- setNames(
  cols$color[cols$Category == "Sex"],
  cols$sample[cols$Category == "Sex"]
)

# Metadata ordered by tree
meta_tree <- meta[match(tree$labels, meta$BioSample), ]

ann_sex   <- unname(sex_col[meta_tree$Sex])
ann_period <- unname(period_col[meta_tree$Period])
ann_tissue <- unname(tissue_col[meta_tree$Tissue])

# Plot
p <- ggtree(tree, branch.length = "none") +
  layout_dendrogram() +
  geom_tippoint(aes(x = x + 0.1),
                color = ann_sex,
                shape = 15,
                size = 0.1) +
  geom_tippoint(aes(x = x + 0.8),
                color = ann_period,
                shape = 15,
                size = 0.1) +
  geom_tippoint(aes(x = x + 1.5),
                color = ann_tissue,
                shape = 15,
                size = 0.1) +
  scale_color_identity() +
  theme(
    legend.position = "none",
    text = element_text(size = 7)
  )

ggsave(
  out_file,
  p,
  width = 15,
  height = 4,
  device = cairo_pdf
)

cat("Samples:", ncol(mat), "\n")
cat("Genes used:", nrow(mat), "\n")
cat("Output:", out_file, "\n")