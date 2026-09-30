#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript period_specific_expression.R <TPM> <counts> <metadata> <tissues.txt> <output_dir>")
}

tpm_file <- args[1]
counts_file <- args[2]
meta_file <- args[3]
tissues_file <- args[4]
outdir <- args[5]

# Parameters from the original workflow
count_threshold <- 6
tpm_threshold <- 0.1
sample_frac_threshold <- 0.20
elbow_threshold <- 0.001
max_tech_pcs <- 15
period_fc_threshold <- 2
period_fdr_threshold <- 0.05

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "period_results"), showWarnings = FALSE)
dir.create(file.path(outdir, "corrected_period_data"), showWarnings = FALSE)

clean_names <- function(x) {
  x <- trimws(x)
  x <- gsub('^"|"$', '', x)
  gsub("^'|'$", '', x)
}

select_pcs_by_elbow <- function(var_exp, d = elbow_threshold) {
  if (length(var_exp) < 3) return(seq_len(min(3, length(var_exp))))
  selected <- integer(0)
  for (n in seq_len(length(var_exp) - 2)) {
    if ((var_exp[n] - var_exp[n + 1]) > d ||
        (var_exp[n] - var_exp[n + 2]) > 2 * d) {
      selected <- c(selected, n)
    } else break
  }
  if (length(selected) == 0) selected <- seq_len(min(3, length(var_exp)))
  head(selected, max_tech_pcs)
}

# Read and filter expression data
tpm_dt <- fread(tpm_file, data.table = FALSE, check.names = FALSE)
count_dt <- fread(counts_file, data.table = FALSE, check.names = FALSE)
meta <- fread(meta_file, data.table = FALSE)

if (!identical(tpm_dt[[1]], count_dt[[1]])) stop("TPM and count gene IDs/order differ.")
if (!all(c("BioSample", "Tissue", "Sex", "Period") %in% colnames(meta))) {
  stop("Metadata must contain BioSample, Tissue, Sex and Period.")
}

genes <- tpm_dt[[1]]
tpm <- as.matrix(tpm_dt[, -1, drop = FALSE])
counts <- as.matrix(count_dt[, -1, drop = FALSE])
storage.mode(tpm) <- "numeric"
storage.mode(counts) <- "numeric"
colnames(tpm) <- clean_names(colnames(tpm))
colnames(counts) <- clean_names(colnames(counts))
meta$BioSample <- clean_names(meta$BioSample)

common <- Reduce(intersect, list(colnames(tpm), colnames(counts), meta$BioSample))
if (length(common) == 0) stop("No matched samples.")
tpm <- tpm[, common, drop = FALSE]
counts <- counts[, common, drop = FALSE]
meta <- meta[match(common, meta$BioSample), , drop = FALSE]

min_samples <- ceiling(sample_frac_threshold * length(common))
keep <- rowSums(tpm >= tpm_threshold, na.rm = TRUE) >= min_samples &
        rowSums(counts >= count_threshold, na.rm = TRUE) >= min_samples
expr_filtered <- tpm[keep, , drop = FALSE]
rownames(expr_filtered) <- genes[keep]

avail_tissues <- if (file.exists(tissues_file)) {
  read.table(tissues_file, header = FALSE, stringsAsFactors = FALSE)$V1
} else {
  names(which(table(meta$Tissue) >= 10))
}
tissues <- intersect(names(which(table(meta$Tissue) >= 10)), avail_tissues)

all_summary <- list()
sum_i <- 1

for (tissue in tissues) {
  cat("Processing tissue:", tissue, "\n")

  tissue_samples <- intersect(meta$BioSample[meta$Tissue == tissue], colnames(expr_filtered))
  if (length(tissue_samples) < 10) next

  meta_t <- meta[match(tissue_samples, meta$BioSample), , drop = FALSE]
  expr_t <- log2(expr_filtered[, tissue_samples, drop = FALSE] + 0.25)
  periods <- unique(meta_t$Period[!is.na(meta_t$Period)])
  if (length(periods) < 2) next

  corrected <- list()
  sample_counts <- list()

  for (period in periods) {
    ps <- meta_t$BioSample[meta_t$Period == period]
    ps <- intersect(ps, colnames(expr_t))
    if (length(ps) < 5) next

    x <- expr_t[, ps, drop = FALSE]
    m <- meta_t[match(ps, meta_t$BioSample), , drop = FALSE]
    x_corrected <- x

    if (length(ps) >= 10) {
      gene_var <- apply(x, 1, var, na.rm = TRUE)
      x_pca <- x[gene_var >= 1e-6, , drop = FALSE]

      if (nrow(x_pca) > 1) {
        pca <- tryCatch(prcomp(t(x_pca), center = TRUE, scale. = TRUE), error = function(e) NULL)
        if (!is.null(pca)) {
          var_exp <- pca$sdev^2 / sum(pca$sdev^2)
          pcs <- select_pcs_by_elbow(var_exp)
          covars <- pca$x[, pcs, drop = FALSE]
          batch_sex <- if (length(unique(m$Sex)) > 1) factor(m$Sex) else NULL
          x_corrected <- removeBatchEffect(x_pca, batch = batch_sex, covariates = covars)
        }
      }
    }

    corrected[[as.character(period)]] <- x_corrected
    sample_counts[[as.character(period)]] <- length(ps)
  }

  if (length(corrected) < 2) next

  tissue_dir <- file.path(outdir, "period_results", tissue)
  dir.create(tissue_dir, recursive = TRUE, showWarnings = FALSE)
  period_results <- list()

  for (target_period in names(corrected)) {
    others <- setdiff(names(corrected), target_period)
    all_genes <- Reduce(intersect, lapply(corrected[c(target_period, others)], rownames))
    if (length(all_genes) < 1000) next

    target_data <- corrected[[target_period]][all_genes, , drop = FALSE]
    other_data <- do.call(cbind, lapply(corrected[others], function(x) x[all_genes, , drop = FALSE]))
    if (ncol(other_data) < 10) next

    expr_combined <- cbind(target_data, other_data)
    group <- factor(c(rep("target", ncol(target_data)), rep("other", ncol(other_data))),
                    levels = c("other", "target"))
    design <- model.matrix(~ 0 + group)
    colnames(design) <- c("other", "target")
    contrast <- makeContrasts(target - other, levels = design)

    fit <- eBayes(contrasts.fit(lmFit(expr_combined, design), contrast))
    result <- topTable(fit, coef = 1, number = Inf, adjust.method = "fdr")
    result$gene_id <- rownames(result)
    result$is_period_specific <- result$adj.P.Val < period_fdr_threshold &
                                 result$logFC > period_fc_threshold
    result <- result[, c("gene_id", setdiff(colnames(result), "gene_id"))]

    fwrite(result,
           file.path(tissue_dir, paste0(tissue, "_", target_period, "_period_specific.tsv")),
           sep = "\t")
    period_results[[target_period]] <- result

    all_summary[[sum_i]] <- data.frame(
      Tissue = tissue,
      Period = target_period,
      Target_Samples = ncol(target_data),
      Other_Samples = ncol(other_data),
      Total_Genes = nrow(result),
      Specific_Genes = sum(result$is_period_specific, na.rm = TRUE)
    )
    sum_i <- sum_i + 1
  }

  saveRDS(
    list(tissue = tissue, period_corrected_data = corrected, meta = meta_t),
    file.path(outdir, "corrected_period_data", paste0(tissue, "_corrected.rds"))
  )
}

if (length(all_summary) > 0) {
  fwrite(rbindlist(all_summary), file.path(outdir, "period_analysis_statistics.tsv"), sep = "\t")
}

cat("Period-specific expression analysis finished.\n")
