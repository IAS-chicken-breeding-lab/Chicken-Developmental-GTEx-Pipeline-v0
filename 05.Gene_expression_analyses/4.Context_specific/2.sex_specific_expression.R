#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(SmartSVA)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript sex_specific_expression.R <TPM> <counts> <gene_annotation> <metadata> <output_dir>")
}

tpm_file   <- args[1]
count_file <- args[2]
gene_file  <- args[3]
meta_file  <- args[4]
outdir     <- args[5]

# Parameters used in the original workflow
fdr_cutoff <- 0.05
large_effect_cutoff <- 1
tpm_threshold <- 0.1
count_threshold <- 6
sample_frac_threshold <- 0.20
target_chromosomes <- c(as.character(1:39), "Z")
set.seed(1)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "Original_TMM"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "Batch_effect_adjusted_TMM"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "ALL"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "FDR"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "Bonferroni"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "Large_effect_FDR"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "Large_effect_Bonferroni"), recursive = TRUE, showWarnings = FALSE)

# Read data
tpm_dt <- fread(tpm_file, data.table = FALSE, check.names = FALSE)
count_dt <- fread(count_file, data.table = FALSE, check.names = FALSE)
gene_info <- fread(gene_file, data.table = FALSE)
meta <- fread(meta_file, data.table = FALSE)

if (!identical(tpm_dt[[1]], count_dt[[1]])) {
  stop("Gene order in TPM and count matrices does not match.")
}

required_meta <- c("BioSample", "Tissue", "Sex", "Period")
if (!all(required_meta %in% colnames(meta))) {
  stop("Missing metadata columns: ", paste(setdiff(required_meta, colnames(meta)), collapse = ", "))
}

if (!all(c("chr", "gene_id") %in% colnames(gene_info))) {
  stop("Gene annotation must contain columns: gene_id and chr")
}

genes <- tpm_dt[[1]]
tpm <- as.matrix(tpm_dt[, -1, drop = FALSE])
counts <- as.matrix(count_dt[, -1, drop = FALSE])
storage.mode(tpm) <- "numeric"
storage.mode(counts) <- "numeric"

# Match samples
common <- Reduce(intersect, list(colnames(tpm), colnames(counts), meta$BioSample))
if (length(common) == 0) stop("No matched samples between expression matrices and metadata.")

tpm <- tpm[, common, drop = FALSE]
counts <- counts[, common, drop = FALSE]
meta <- meta[match(common, meta$BioSample), , drop = FALSE]

# Global TPM/count/chromosome filtering
min_samples <- ceiling(sample_frac_threshold * length(common))
keep <- rowSums(tpm > tpm_threshold, na.rm = TRUE) >= min_samples &
        rowSums(counts >= count_threshold, na.rm = TRUE) >= min_samples &
        genes %in% gene_info$gene_id[gene_info$chr %in% target_chromosomes]

genes <- genes[keep]
counts <- counts[keep, , drop = FALSE]
rownames(counts) <- genes

cat("Samples:", ncol(counts), "\n")
cat("Genes after filtering:", nrow(counts), "\n")

summary_list <- list()
idx_summary <- 1

for (period in unique(meta$Period)) {
  for (tissue in unique(meta$Tissue)) {

    sample_idx <- which(meta$Period == period & meta$Tissue == tissue)
    if (length(sample_idx) == 0) next

    meta_sub <- meta[sample_idx, , drop = FALSE]
    count_sub <- counts[, meta_sub$BioSample, drop = FALSE]

    meta_sub$Sex <- factor(meta_sub$Sex, levels = c("F", "M"))
    meta_sub <- meta_sub[!is.na(meta_sub$Sex), , drop = FALSE]
    count_sub <- count_sub[, meta_sub$BioSample, drop = FALSE]

    if (nlevels(droplevels(meta_sub$Sex)) < 2) next

    cat("Processing", period, tissue, "n =", ncol(count_sub), "\n")

    # edgeR filtering + TMM normalization
    y <- DGEList(counts = count_sub, group = meta_sub$Sex)
    y <- y[filterByExpr(y), , keep.lib.sizes = FALSE]
    if (nrow(y) == 0) next
    y <- calcNormFactors(y, method = "TMM")
    count_norm <- as.data.frame(cpm(y))

    # SmartSVA
    mod <- model.matrix(~ meta_sub$Sex)
    rownames(mod) <- colnames(count_norm)
    n_sv <- num.sv(count_norm, mod)

    if (n_sv > 0) {
      sv <- smartsva.cpp(as.matrix(count_norm), mod = mod, n.sv = n_sv)$sv
      rownames(sv) <- colnames(count_norm)
      count_adjusted <- removeBatchEffect(as.matrix(count_norm), covariates = sv)
    } else {
      count_adjusted <- as.matrix(count_norm)
    }

    prefix <- paste(period, tissue, sep = "_")
    write.table(count_norm,
                file.path(outdir, "Original_TMM", paste0(prefix, "_TMM.tsv")),
                sep = "\t", quote = FALSE)
    write.table(count_adjusted,
                file.path(outdir, "Batch_effect_adjusted_TMM", paste0(prefix, "_adjusted_expression.tsv")),
                sep = "\t", quote = FALSE)

    # Wilcoxon test and Male/Female fold change
    pvalues <- apply(count_adjusted, 1, function(x) {
      wilcox.test(x[meta_sub$Sex == "M"], x[meta_sub$Sex == "F"])$p.value
    })

    mean_f <- rowMeans(count_adjusted[, meta_sub$Sex == "F", drop = FALSE])
    mean_m <- rowMeans(count_adjusted[, meta_sub$Sex == "M", drop = FALSE])
    log2fc <- log2(mean_m / mean_f)

    result <- data.frame(
      Gene = rownames(count_adjusted),
      log2foldChange = log2fc,
      pValue = pvalues,
      FDR = p.adjust(pvalues, method = "fdr"),
      Bonferroni = p.adjust(pvalues, method = "bonferroni"),
      Chromosome = gene_info$chr[match(rownames(count_adjusted), gene_info$gene_id)],
      Direction = ifelse(log2fc > 0, "Male biased", "Female biased"),
      stringsAsFactors = FALSE
    )
    result <- na.omit(result)

    fwrite(result, file.path(outdir, "ALL", paste0("ALL_", prefix, ".tsv")), sep = "\t")
    fwrite(result[result$FDR < fdr_cutoff, ],
           file.path(outdir, "FDR", paste0("FDR_", prefix, ".tsv")), sep = "\t")
    fwrite(result[result$Bonferroni < 0.05, ],
           file.path(outdir, "Bonferroni", paste0("Bonferroni_", prefix, ".tsv")), sep = "\t")
    fwrite(result[result$FDR < fdr_cutoff & abs(result$log2foldChange) > large_effect_cutoff, ],
           file.path(outdir, "Large_effect_FDR", paste0("Large_effect_FDR_", prefix, ".tsv")), sep = "\t")
    fwrite(result[result$Bonferroni < 0.05 & abs(result$log2foldChange) > large_effect_cutoff, ],
           file.path(outdir, "Large_effect_Bonferroni", paste0("Large_effect_Bonferroni_", prefix, ".tsv")), sep = "\t")

    summary_list[[idx_summary]] <- data.frame(
      Period = period,
      Tissue = tissue,
      Samples = ncol(count_sub),
      Tested_Genes = nrow(result),
      SVs = n_sv,
      FDR_Significant = sum(result$FDR < fdr_cutoff),
      Bonferroni_Significant = sum(result$Bonferroni < 0.05)
    )
    idx_summary <- idx_summary + 1
  }
}

if (length(summary_list) > 0) {
  fwrite(rbindlist(summary_list), file.path(outdir, "summary_statistics.tsv"), sep = "\t")
}

cat("Sex-specific expression analysis finished.\n")
