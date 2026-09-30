#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(variancePartition)
})

args <- commandArgs(trailingOnly = TRUE)

expr_file <- args[1]
meta_file <- args[2]
out_file  <- args[3]

# -------------------------
# Read data
# -------------------------
expr <- fread(expr_file, data.table = FALSE)
rownames(expr) <- expr[[1]]
expr <- expr[, -1, drop = FALSE]

meta <- fread(meta_file, data.table = FALSE)
sample_id_col <- colnames(meta)[1]

required_cols <- c(
  "Individual", "Tissue", "Sex", "Stage",
  "Weight", "GHRH", "INS", "GH",
  "CORT", "IGF-1", "Ghrelin"
)

if (!all(required_cols %in% colnames(meta))) {
  stop("Missing required columns: ",
       paste(setdiff(required_cols, colnames(meta)), collapse = ", "))
}

# -------------------------
# Expression filtering
# -------------------------
expr <- expr[rowSums(expr >= 0.1) > 0, , drop = FALSE]

keep_samples <- colMeans(expr >= 0.1) >= 0.20
expr <- expr[, keep_samples, drop = FALSE]

expr <- expr[rowSums(expr > 0) > 0, , drop = FALSE]

# -------------------------
# Match samples
# -------------------------
common_samples <- intersect(colnames(expr), meta[[sample_id_col]])

if (length(common_samples) == 0)
  stop("No matched samples between expression and metadata.")

expr <- expr[, common_samples, drop = FALSE]
meta <- meta[match(common_samples, meta[[sample_id_col]]), , drop = FALSE]

# -------------------------
# Remove missing covariates
# -------------------------
model_vars <- c(
  "Individual", "Tissue", "Sex", "Stage",
  "Weight", "GHRH", "INS", "GH",
  "cortisol", "IGF", "ghrelin"
)

keep <- complete.cases(meta[, model_vars])
meta <- meta[keep, , drop = FALSE]
expr <- expr[, keep, drop = FALSE]

# -------------------------
# Transform variables
# -------------------------
continuous_vars <- c(
  "Weight", "GHRH", "INS", "GH",
  "cortisol", "IGF", "ghrelin"
)

meta[, continuous_vars] <- scale(meta[, continuous_vars])

meta$Individual <- factor(meta$Individual)
meta$Tissue <- factor(meta$Tissue)
meta$Sex <- factor(meta$Sex)
meta$Stage <- factor(meta$Stage)

expr <- log2(as.matrix(expr) + 1)

# -------------------------
# Variance partitioning
# -------------------------
form <- ~ Weight + GHRH + INS + GH + cortisol + IGF + ghrelin +
  (1 | Tissue) +
  (1 | Sex) +
  (1 | Stage) +
  (1 | Individual)

varPart <- fitExtractVarPartModel(
  expr,
  form,
  meta
)

write.csv(varPart, out_file, row.names = TRUE)

cat("Samples:", ncol(expr), "\n")
cat("Genes:", nrow(expr), "\n")
cat("Results:", out_file, "\n")