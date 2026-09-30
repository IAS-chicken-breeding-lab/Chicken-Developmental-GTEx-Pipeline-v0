#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript tissue_specific_expression.R <TPM> <counts> <metadata> <tissues.txt> <output_dir>")
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
tissue_fc_threshold <- 2
tissue_fdr_threshold <- 0.05

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "corrected_data"), showWarnings = FALSE)
dir.create(file.path(outdir, "tissue_results"), showWarnings = FALSE)

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

safe_prcomp <- function(x) {
  sds <- apply(x, 2, sd, na.rm = TRUE)
  x <- x[, sds > 1e-8, drop = FALSE]
  if (ncol(x) < 3) return(NULL)
  tryCatch(prcomp(x, center = TRUE, scale. = TRUE),
           error = function(e) tryCatch(prcomp(x, center = TRUE, scale. = FALSE),
                                        error = function(e2) NULL))
}

# Read and filter data
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

# Step 1: correct expression within each tissue
all_corrected <- list()

for (tissue in tissues) {
  cat("Correcting tissue:", tissue, "\n")

  samples <- intersect(meta$BioSample[meta$Tissue == tissue], colnames(expr_filtered))
  if (length(samples) < 5) next

  meta_t <- meta[match(samples, meta$BioSample), , drop = FALSE]
  x <- log2(expr_filtered[, samples, drop = FALSE] + 0.25)

  pca_covars <- NULL
  if (length(samples) >= 10) {
    row_var <- apply(x, 1, var, na.rm = TRUE)
    x_pca <- x[row_var > 1e-10, , drop = FALSE]
    if (nrow(x_pca) >= 3) {
      pca <- safe_prcomp(t(x_pca))
      if (!is.null(pca)) {
        var_exp <- pca$sdev^2 / sum(pca$sdev^2)
        pcs <- select_pcs_by_elbow(var_exp)
        pca_covars <- pca$x[, pcs, drop = FALSE]
      }
    }
  }

  batch_period <- if (length(unique(meta_t$Period)) > 1) factor(meta_t$Period) else NULL
  batch_sex <- if (length(unique(meta_t$Sex)) > 1) factor(meta_t$Sex) else NULL

  x_corrected <- x
  if (!is.null(batch_period) || !is.null(batch_sex) || !is.null(pca_covars)) {
    x_corrected <- tryCatch(
      removeBatchEffect(x, batch = batch_period, batch2 = batch_sex, covariates = pca_covars),
      error = function(e) x
    )
  }

  all_corrected[[tissue]] <- list(
    tissue = tissue,
    samples = samples,
    expr_corrected = x_corrected,
    meta = meta_t,
    pca_covars_used = ifelse(is.null(pca_covars), 0, ncol(pca_covars))
  )

  saveRDS(all_corrected[[tissue]],
          file.path(outdir, "corrected_data", paste0(tissue, ".rds")))
}

# Step 2: target tissue vs all other tissues
summary_list <- list()
sum_i <- 1

for (tissue in names(all_corrected)) {
  cat("Testing tissue:", tissue, "\n")

  target <- all_corrected[[tissue]]$expr_corrected
  others <- all_corrected[setdiff(names(all_corrected), tissue)]
  if (length(others) == 0) next

  common_genes <- Reduce(intersect, c(list(rownames(target)), lapply(others, function(x) rownames(x$expr_corrected))))
  if (length(common_genes) == 0) next

  target <- target[common_genes, , drop = FALSE]
  other_data <- do.call(cbind, lapply(others, function(x) x$expr_corrected[common_genes, , drop = FALSE]))
  if (ncol(other_data) < 20) next

  expr_combined <- cbind(target, other_data)
  group <- factor(c(rep("target", ncol(target)), rep("other", ncol(other_data))),
                  levels = c("other", "target"))
  design <- model.matrix(~ 0 + group)
  colnames(design) <- c("other", "target")
  contrast <- makeContrasts(target - other, levels = design)

  fit <- eBayes(contrasts.fit(lmFit(expr_combined, design), contrast))
  result <- topTable(fit, coef = 1, number = Inf, adjust.method = "fdr")
  result$gene_id <- rownames(result)
  result$is_tissue_specific <- result$adj.P.Val < tissue_fdr_threshold &
                               result$logFC > tissue_fc_threshold
  result <- result[, c("gene_id", setdiff(colnames(result), "gene_id"))]

  fwrite(result,
         file.path(outdir, "tissue_results", paste0(tissue, "_tissue_specific.tsv")),
         sep = "\t")

  summary_list[[sum_i]] <- data.frame(
    Tissue = tissue,
    Target_Samples = ncol(target),
    Other_Samples = ncol(other_data),
    Total_Genes = nrow(result),
    Specific_Genes = sum(result$is_tissue_specific, na.rm = TRUE),
    PCA_Covars_Used = all_corrected[[tissue]]$pca_covars_used
  )
  sum_i <- sum_i + 1
}

if (length(summary_list) > 0) {
  fwrite(rbindlist(summary_list), file.path(outdir, "tissue_analysis_statistics.tsv"), sep = "\t")
}

cat("Tissue-specific expression analysis finished.\n")
