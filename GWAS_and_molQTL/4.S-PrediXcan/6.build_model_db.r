#!/usr/bin/env Rscript
#
# build_model_db.r
# Build SQLite prediction model databases (PredictDB format) from
# combined TWAS model outputs.
#
# Usage:
#   Rscript build_model_db.r \
#     --combined-dir <dir> \
#     --out-dir <dir> \
#     [--tissue-pattern <regex>] \
#     [--chromosomes <list>] \
#     [--zscore-pval <float>] \
#     [--rho-cutoff <float>]
#
# Arguments:
#   --combined-dir    Directory containing CombinedModel/{summary,weights}/<tissue>/
#   --out-dir         Directory for output .db files
#
# Options:
#   --tissue-pattern  Regex to select tissues (default: "^A")
#   --chromosomes     Comma-separated chromosome list (default: 1..39)
#   --zscore-pval     z-score p-value threshold for filtered DB (default: 0.05)
#   --rho-cutoff      rho average cutoff for filtered DB (default: 0.1)
#

suppressPackageStartupMessages({
  library(RSQLite)
  library(dplyr)
})

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)

get_opt <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i + 1 > length(args)) stop("Missing value for ", flag)
  args[i + 1]
}

combined_dir   <- get_opt("--combined-dir")
out_dir        <- get_opt("--out-dir")
tissue_pattern <- get_opt("--tissue-pattern", "^A")
chromosomes    <- get_opt("--chromosomes", NULL)
zscore_pval    <- as.numeric(get_opt("--zscore-pval", 0.05))
rho_cutoff     <- as.numeric(get_opt("--rho-cutoff", 0.1))

if (is.null(combined_dir) || is.null(out_dir)) {
  stop("Required options: --combined-dir and --out-dir")
}

if (!dir.exists(combined_dir)) {
  stop("Combined directory not found: ", combined_dir)
}

summary_dir <- file.path(combined_dir, "summary")
weights_dir <- file.path(combined_dir, "weights")

if (!dir.exists(summary_dir)) {
  stop("Summary directory not found: ", summary_dir)
}
if (!dir.exists(weights_dir)) {
  stop("Weights directory not found: ", weights_dir)
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Chromosome list
if (!is.null(chromosomes)) {
  chr_list <- as.integer(unlist(strsplit(chromosomes, ",")))
} else {
  chr_list <- 1:39
}

# Start
cat("Building chicken eQTL model databases\n")
cat("Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Summary dir     :", summary_dir, "\n")
cat("Weights dir     :", weights_dir, "\n")
cat("Output dir      :", out_dir, "\n")
cat("Tissue pattern  :", tissue_pattern, "\n")
cat("Chromosomes     :", paste(range(chr_list), collapse = "-"), "\n")

tissues <- list.files(summary_dir, pattern = tissue_pattern)
if (length(tissues) == 0) {
  stop("No tissues matching pattern '", tissue_pattern, "' found.")
}

cat("Tissues found:\n  ", paste(tissues, collapse = ", "), "\n\n")

# Helper: rename columns safely
safe_rename <- function(df, map) {
  for (from in names(map)) {
    to <- map[[from]]
    if (from %in% names(df)) {
      df <- df %>% rename(!!to := !!sym(from))
    } else {
      df[[to]] <- NA
    }
  }
  df
}

# Helper: load summaries or weights for a tissue
load_tissue_files <- function(dir, tissue, suffix, chr_list) {
  all <- list()
  for (chr in chr_list) {
    f <- file.path(dir, tissue,
                   sprintf("%s_Model_training_chr%d_%s.txt", tissue, chr, suffix))
    if (file.exists(f)) {
      all[[length(all) + 1]] <- read.table(f, header = TRUE, sep = "\t",
                                           stringsAsFactors = FALSE)
    }
  }
  if (length(all) == 0) return(NULL)
  bind_rows(all)
}

# Main loop

for (tissue in tissues) {
  cat("Tissue:", tissue, "\n")

  #  Read model summaries 
  model_all <- load_tissue_files(summary_dir, tissue, "model_summaries", chr_list)
  if (is.null(model_all) || nrow(model_all) == 0) {
    cat("  [SKIP] No model_summaries data\n")
    next
  }
  cat("  model_summaries rows:", nrow(model_all), "\n")

  if ("gene_id" %in% names(model_all)) {
    model_all <- rename(model_all, gene = gene_id)
  }

  #  Read weights 
  weights_all <- load_tissue_files(weights_dir, tissue, "weights", chr_list)
  if (is.null(weights_all) || nrow(weights_all) == 0) {
    cat("  [SKIP] No weights data\n")
    next
  }
  cat("  weights rows:", nrow(weights_all), "\n")

  if ("gene_id" %in% names(weights_all)) {
    weights_all <- rename(weights_all, gene = gene_id)
  }

  #  Full database 
  full_db <- file.path(out_dir,
                       sprintf("Chicken_%s_ElasticNet_models.db", tissue))
  full_conn <- dbConnect(SQLite(), full_db)
  dbWriteTable(full_conn, "model_summaries", model_all, overwrite = TRUE)
  dbWriteTable(full_conn, "weights", weights_all, overwrite = TRUE)
  dbWriteTable(full_conn, "meta_data",
               data.frame(key = c("db_type", "population", "tissue"),
                          value = c("PredictionModel", "Chicken", tissue)),
               overwrite = TRUE)
  dbDisconnect(full_conn)
  cat("  Full DB:", basename(full_db),
      sprintf("(%.1f MB)", file.size(full_db) / 1024 / 1024), "\n")

  #  Filtered database 
  in_conn <- dbConnect(SQLite(), full_db)
  filtered_db <- file.path(out_dir,
                           sprintf("Chicken_%s_ElasticNet_models_filtered_signif.db", tissue))
  out_conn <- dbConnect(SQLite(), filtered_db)

  q <- sprintf("SELECT * FROM model_summaries WHERE zscore_pval < %s AND rho_avg > %s",
               zscore_pval, rho_cutoff)
  msf <- dbGetQuery(in_conn, q)
  cat("  Significant models:", nrow(msf), "\n")

  if (nrow(msf) == 0) {
    dbExecute(out_conn,
              "CREATE TABLE extra(gene TEXT, genename TEXT, n.snps.in.model INTEGER)")
    dbExecute(out_conn,
              "CREATE TABLE weights(gene TEXT, rsid TEXT, weight REAL, eff_allele TEXT, ref_allele TEXT)")
    cat("  [WARN] No significant models; empty structure created\n")
  } else {
    rename_map <- c(
      "gene_name"             = "genename",
      "n_snps_in_model"       = "n.snps.in.model",
      "n_snps_in_window"      = "n.snps.in.window",
      "in_sample_R2"          = "pred.perf.R2",
      "nested_cv_fisher_pval" = "pred.perf.pval",
      "test_R2_avg"           = "test.r2.avg",
      "test_R2_sd"            = "test.r2.sd",
      "cv_R2_avg"             = "cv.r2.avg",
      "cv_R2_sd"              = "cv.r2.sd",
      "rho_avg"               = "rho.avg",
      "rho_se"                = "rho.se",
      "rho_zscore"            = "rho.zscore",
      "zscore_estimate"       = "cv.zscore.est",
      "zscore_pval"           = "cv.pval.est"
    )
    extra <- safe_rename(msf, rename_map)
    extra <- extra %>% mutate(pred.perf.qval = NA)
    dbWriteTable(out_conn, "extra", extra, overwrite = TRUE)

    wq <- sprintf("SELECT * FROM weights WHERE gene IN (
                     SELECT gene FROM model_summaries
                     WHERE zscore_pval < %s AND rho_avg > %s)",
                  zscore_pval, rho_cutoff)
    wf <- dbGetQuery(in_conn, wq)

    if ("beta" %in% names(wf)) {
      wf <- wf %>% rename(weight = beta)
    }
    if (!"eff_allele" %in% names(wf)) {
      wf$eff_allele <- if ("alt" %in% names(wf)) wf$alt else NA
    }
    if (!"ref_allele" %in% names(wf)) {
      wf$ref_allele <- if ("ref" %in% names(wf)) wf$ref else NA
    }
    wf <- wf %>% select(gene, rsid, weight, eff_allele, ref_allele)

    dbWriteTable(out_conn, "weights", wf, overwrite = TRUE)
    cat("  Filtered DB written\n")
  }

  dbWriteTable(out_conn, "meta_data",
               data.frame(key = c("db_type", "population", "tissue"),
                          value = c("PredictionModel", "Chicken", tissue)),
               overwrite = TRUE)
  dbDisconnect(in_conn)
  dbDisconnect(out_conn)
  cat("  Filtered DB:", basename(filtered_db),
      sprintf("(%.1f MB)", file.size(filtered_db) / 1024 / 1024), "\n")
}

cat("All tissues processed. Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")