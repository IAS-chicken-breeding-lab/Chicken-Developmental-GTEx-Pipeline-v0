#!/usr/bin/env Rscript
#
# parallel_tiss_chrom_training.r
# Train Elastic Net TWAS prediction models for one tissue and one chromosome.
#
# Usage:
#   Rscript parallel_tiss_chrom_training.r <chrom> <tissue> <type> \
#     --base-dir <dir> \
#     --covar-dir <dir> \
#     --gene-annot <file> \
#     [--nested-cv-script <path>] \
#     [--max-cores <n>]
#
# Arguments:
#   chrom            Chromosome (string, e.g., "1")
#   tissue           Tissue identifier (e.g., "A_AF")
#   type             Data type label (e.g., "eQTL", "sQTL")
#
# Options:
#   --base-dir         Base directory containing model training inputs
#   --covar-dir        Directory with covariate files (<tissue>_cov.txt)
#   --gene-annot       Gene annotation file
#   --nested-cv-script R script implementing the nested CV Elastic Net model
#                      (default: ./parallel_nested_cv_elnet.r)
#   --max-cores        Maximum number of parallel cores (default: 10)
#

suppressPackageStartupMessages({
  library(parallel)
  library(dplyr)
  library(glmnet)
  library(reshape2)
  library(methods)
  library(data.table)
})

`%&%` <- function(a, b) paste(a, b, sep = "")

# Parse command-line arguments
argv <- commandArgs(trailingOnly = TRUE)

if (length(argv) < 3) {
  cat("Usage: Rscript parallel_tiss_chrom_training.r <chrom> <tissue> <type>",
      "--base-dir <dir> --covar-dir <dir> --gene-annot <file>",
      "[--nested-cv-script <path>] [--max-cores <n>]\n")
  quit(status = 1)
}

chrom  <- argv[1]
tissue <- argv[2]
type   <- argv[3]

get_opt <- function(flag, default = NULL) {
  i <- which(argv == flag)
  if (length(i) == 0) return(default)
  if (i + 1 > length(argv)) stop("Missing value for ", flag)
  argv[i + 1]
}

base_dir         <- get_opt("--base-dir")
covar_dir        <- get_opt("--covar-dir")
gene_annot_file  <- get_opt("--gene-annot")
nested_cv_script <- get_opt("--nested-cv-script", "./parallel_nested_cv_elnet.r")
max_cores        <- as.integer(get_opt("--max-cores", 10))

if (is.null(base_dir) || is.null(covar_dir) || is.null(gene_annot_file)) {
  stop("Missing required options: --base-dir, --covar-dir, --gene-annot")
}

if (!file.exists(nested_cv_script)) {
  stop("Nested CV script not found: ", nested_cv_script)
}

# Derived paths
expression_file <- file.path(base_dir, "tissue_expression",
                             paste0(tissue, ".", type, ".transformed_expression.txt"))
covariates_file <- file.path(covar_dir, paste0(tissue, "_cov.txt"))
snp_annot_file  <- file.path(base_dir, "tissue_snp_annotation", tissue,
                             paste0(tissue, ".snp_annot.chr", chrom, ".txt"))
genotype_file   <- file.path(base_dir, "tissue_gene_matrix", tissue,
                             paste0(tissue, ".genotype.chr", chrom, ".txt"))
prefix <- "Model_training"

# Logging
cat("=== TWAS model training ===\n")
cat("Tissue:", tissue, " Chromosome:", chrom, " Type:", type, "\n\n")

check_file <- function(f, desc) {
  if (!file.exists(f)) stop(desc, " not found: ", f)
  cat(desc, ":", f, "\n")
}

check_file(gene_annot_file, "Gene annotation file")
check_file(expression_file, "Expression file")
check_file(covariates_file, "Covariates file")
check_file(snp_annot_file, "SNP annotation file")
check_file(genotype_file, "Genotype file")

# Helper: read gene annotation for a chromosome
get_gene_annotation <- function(gene_annot_file_name, chrom) {
  gene_df <- fread(gene_annot_file_name, header = TRUE, sep = "\t")
  gene_df$chr <- as.character(gene_df$chr)
  chrom_chr <- as.character(chrom)
  gene_df <- gene_df[gene_df$chr == chrom_chr, ]
  cat("Genes on chromosome", chrom_chr, ":", nrow(gene_df), "\n")
  if (!"gene_name" %in% colnames(gene_df)) {
    stop("Gene annotation file must contain a 'gene_name' column.")
  }
  gene_df
}

# Helper: read expression matrix
get_gene_expression <- function(gene_expression_file_name, gene_annot) {
  cat("Reading expression matrix ...\n")
  expr_df <- fread(gene_expression_file_name, header = TRUE, sep = "\t",
                   check.names = FALSE)
  expr_df <- as.data.frame(expr_df)

  cat("Expression matrix dimensions:", nrow(expr_df), "x", ncol(expr_df), "\n")
  cat("First column:", names(expr_df)[1], "\n")

  # Sample IDs
  id_vec <- as.character(expr_df[[1]])
  if (any(id_vec == "" | is.na(id_vec) | duplicated(id_vec))) {
    cat("Warning: duplicated or empty sample IDs, applying make.unique ...\n")
    id_vec[id_vec == "" | is.na(id_vec)] <- "sample"
    id_vec <- make.unique(id_vec)
  }
  rownames(expr_df) <- id_vec
  expr_df <- expr_df[, -1, drop = FALSE]

  expr_genes  <- colnames(expr_df)
  annot_genes <- gene_annot$gene_name

  cat("Genes in expression matrix:", length(expr_genes), "\n")
  cat("Genes in annotation       :", length(annot_genes), "\n")

  # Exact match
  common_genes <- intersect(expr_genes, annot_genes)
  cat("Exact-matched genes:", length(common_genes), "\n")

  # Try stripping version suffixes
  if (length(common_genes) == 0) {
    cat("No exact match; trying version-stripped matching ...\n")
    expr_prefix  <- gsub("\\..*", "", expr_genes)
    annot_prefix <- gsub("\\..*", "", annot_genes)
    common_prefix <- intersect(expr_prefix, annot_prefix)
    cat("Version-stripped matches:", length(common_prefix), "\n")
    if (length(common_prefix) > 0) {
      keep_idx <- which(expr_prefix %in% common_prefix)
      expr_df  <- expr_df[, keep_idx, drop = FALSE]
      colnames(expr_df) <- expr_prefix[keep_idx]
      common_genes <- colnames(expr_df)
    }
  }

  if (length(common_genes) == 0) {
    stop("No matching genes. Check that expression column names match annotation 'gene_name'.")
  }

  expr_final <- expr_df[, common_genes, drop = FALSE]
  cat("Final expression matrix:", nrow(expr_final), "samples x",
      ncol(expr_final), "genes\n")
  expr_final
}

# Helper: read covariates
get_covariates <- function(covariates_file_name) {
  cat("Reading covariates ...\n")
  covar_df <- fread(covariates_file_name, header = TRUE, sep = "\t",
                    check.names = FALSE)
  covar_df <- as.data.frame(covar_df)

  cat("Covariate file dimensions:", nrow(covar_df), "x", ncol(covar_df), "\n")
  cat("First column (should be CovarName):", names(covar_df)[1], "\n")

  covar_names <- as.character(covar_df[[1]])
  if (any(is.na(covar_names) | covar_names == "")) {
    cat("Warning: empty covariate names, auto-generating ...\n")
    empty_idx <- which(is.na(covar_names) | covar_names == "")
    covar_names[empty_idx] <- paste0("Covar", empty_idx)
  }

  covar_t <- as.data.frame(t(covar_df[, -1, drop = FALSE]))
  colnames(covar_t) <- make.names(covar_names, unique = TRUE)

  covar_t$Id <- rownames(covar_t)
  covar_t <- covar_t[, c("Id", setdiff(colnames(covar_t), "Id"))]

  cat("Covariates loaded:", nrow(covar_t), "samples,",
      ncol(covar_t) - 1, "covariates\n")
  cat("Covariate names:", paste(colnames(covar_t)[-1], collapse = ", "), "\n")
  covar_t
}

# Main
tryCatch({
  cat("Loading data ...\n")
  gene_annot <- get_gene_annotation(gene_annot_file, chrom)
  expr_df    <- get_gene_expression(expression_file, gene_annot)
  covariates <- get_covariates(covariates_file)

  n_gene <- ncol(expr_df)
  if (n_gene == 0) stop("No genes found after matching.")

  cat("Processing:", tissue, "chromosome", chrom, "genes:", n_gene, "\n")

  parallel_main <- function(i) {
    source(nested_cv_script)
    main(snp_annot_file, gene_annot, genotype_file, expr_df, covariates,
         as.numeric(chrom), prefix, tissue, type, as.numeric(i),
         null_testing = FALSE)
  }

  num_cores <- max(1, min(max_cores, detectCores() - 1))
  cat("Using", num_cores, "cores for parallel computation\n")
  cl <- makeCluster(num_cores, type = "FORK")

  clusterEvalQ(cl, {
    library(dplyr)
    library(glmnet)
    library(reshape2)
    library(methods)
    library(data.table)
  })

  clusterExport(cl, c("snp_annot_file", "gene_annot", "genotype_file",
                      "expr_df", "covariates", "prefix",
                      "tissue", "type", "chrom", "base_dir",
                      "nested_cv_script"))

  cat("Running parallel analysis ...\n")
  parLapply(cl, seq_len(n_gene), parallel_main)
  stopCluster(cl)

  cat("Done:", tissue, "chromosome", chrom, "\n")
  cat("Results in:", base_dir, "\n")

}, error = function(e) {
  cat("Error in", tissue, "chromosome", chrom, ":", e$message, "\n")
  stop(e)
})

cat("=== TWAS model training finished ===\n")
cat("Time:", date(), "\n")