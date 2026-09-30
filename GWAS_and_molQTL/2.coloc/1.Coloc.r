#!/usr/bin/env Rscript

#
# Usage:
#   Rscript run_coloc.R \
#     --gwas-dir <dir> \
#     --gwas-sample <file> \
#     --eqtl-dir <dir> \
#     --eqtl-sample <file> \
#     --outdir <dir> \
#     [--threads <n>] \
#     [gwas_file1 gwas_file2 ...]

suppressPackageStartupMessages({
  required_packages <- c("data.table", "dplyr", "coloc", "tools",
                         "foreach", "doParallel")
  missing <- required_packages[!sapply(required_packages, requireNamespace,
                                       quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing required packages: ", paste(missing, collapse = ", "))
  }
  library(data.table)
  library(dplyr)
  library(coloc)
  library(tools)
  library(foreach)
  library(doParallel)
})

# Default paths (edit these if you run the script without arguments)
default_gwas_dir      <- "/path/to/gwas_dir"
default_gwas_sample   <- "/path/to/gwas_sample.txt"
default_eqtl_dir      <- "/path/to/eqtl_dir"
default_eqtl_sample   <- "/path/to/eqtl_sample.txt"
default_out_root      <- "/path/to/output_dir"
default_ncores        <- NULL   # NULL = auto-detect

# Helper functions
log_msg <- function(...) {
  msg <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
               " ", paste(..., collapse = " "), "\n")
  message(msg)
  flush.console()
}

safe_name <- function(x) {
  x <- as.character(x)
  x <- gsub("[/\\\\:*?\"<>|]", "_", x)
  x <- gsub("[^A-Za-z0-9._-]", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

parse_eqtl_info <- function(filename) {
  base <- basename(filename)
  parts <- strsplit(base, "_")[[1]]
  if (length(parts) < 2) stop(paste0("Invalid eQTL filename: ", base))
  list(
    Period = toupper(parts[1]),
    Tissue = toupper(parts[2]),
    file   = filename
  )
}

extract_pheno <- function(filename) {
  base <- basename(filename)
  base_noext <- file_path_sans_ext(base)
  pheno <- sub("_extracted(_regions)?$", "", base_noext)
  if (pheno == base_noext) pheno <- base_noext
  pheno
}

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i + 1 > length(args)) stop(paste("Missing value for", flag))
  args[i + 1]
}

gwas_dir      <- get_arg("--gwas-dir",      default_gwas_dir)
gwas_sample_f <- get_arg("--gwas-sample",   default_gwas_sample)
eqtl_dir      <- get_arg("--eqtl-dir",      default_eqtl_dir)
eqtl_sample_f <- get_arg("--eqtl-sample",   default_eqtl_sample)
out_root      <- get_arg("--outdir",        default_out_root)
ncores        <- as.integer(get_arg("--threads", default_ncores))

# Remaining arguments are GWAS filenames
flag_idx <- which(args %in% c("--gwas-dir", "--gwas-sample", "--eqtl-dir",
                              "--eqtl-sample", "--outdir", "--threads"))
to_remove <- sort(unique(c(flag_idx, flag_idx + 1)))
selected_gwas <- if (length(to_remove) > 0) args[-to_remove] else args

if (is.null(gwas_dir) || is.null(gwas_sample_f) ||
    is.null(eqtl_dir) || is.null(eqtl_sample_f) ||
    is.null(out_root)) {
  stop("Missing required paths. Edit the defaults or provide command-line arguments.")
}

dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

log_msg("====== GWAS x eQTL colocalization ======")
log_msg("GWAS dir       : ", gwas_dir)
log_msg("GWAS sample    : ", gwas_sample_f)
log_msg("eQTL dir       : ", eqtl_dir)
log_msg("eQTL sample    : ", eqtl_sample_f)
log_msg("Output dir     : ", out_root)

# Read sample size tables
stopifnot(file.exists(gwas_sample_f))
gwas_sample_df <- fread(gwas_sample_f, header = FALSE,
                        col.names = c("file", "N"))
gwas_sample_df[, file := basename(file)]
log_msg("GWAS sample size records: ", nrow(gwas_sample_df))

stopifnot(file.exists(eqtl_sample_f))
eqtl_sample <- fread(eqtl_sample_f)
required_cols <- c("Period", "Tissue", "SampleSize")
if (!all(required_cols %in% names(eqtl_sample))) {
  stop("eQTL sample file missing columns: ",
       paste(setdiff(required_cols, names(eqtl_sample)), collapse = ", "))
}
eqtl_sample <- eqtl_sample %>%
  mutate(
    Period     = toupper(trimws(Period)),
    Tissue     = toupper(trimws(Tissue)),
    SampleSize = as.integer(SampleSize)
  ) %>%
  filter(!is.na(SampleSize), SampleSize > 0) %>%
  distinct(Period, Tissue, .keep_all = TRUE)
log_msg("eQTL sample size records: ", nrow(eqtl_sample))

# Collect GWAS files
gwas_files_all <- list.files(gwas_dir, pattern = "\\.txt$", full.names = TRUE)
if (length(gwas_files_all) == 0) stop("No GWAS files found in ", gwas_dir)

if (length(selected_gwas) > 0) {
  gwas_files <- gwas_files_all[basename(gwas_files_all) %in% selected_gwas]
  if (length(gwas_files) == 0) {
    stop("None of the specified GWAS files were found.")
  }
  log_msg("Selected GWAS files: ", length(gwas_files))
} else {
  gwas_files <- gwas_files_all
  log_msg("Processing all GWAS files: ", length(gwas_files))
}

gwas_files <- gwas_files[basename(gwas_files) %in% gwas_sample_df$file]
if (length(gwas_files) == 0) stop("No GWAS files with sample sizes.")
log_msg("GWAS files with sample sizes: ", length(gwas_files))

# Collect eQTL files (Period_Tissue_eQTL.txt)
eqtl_files_all <- list.files(eqtl_dir,
                             pattern = "_eQTL\\.txt$",
                             full.names = TRUE)
if (length(eqtl_files_all) == 0) stop("No eQTL files found in ", eqtl_dir)
eqtl_info_list <- lapply(eqtl_files_all, parse_eqtl_info)
log_msg("eQTL files: ", length(eqtl_info_list))

# Parallel setup
if (is.na(ncores) || is.null(ncores)) {
  ncores <- max(1, detectCores() - 2)
}
log_msg("Cores: ", ncores)
cl <- makeCluster(ncores, outfile = "")
registerDoParallel(cl)

# Main loop
foreach(gwas_file = gwas_files,
        .packages = c("data.table", "dplyr", "coloc", "tools"),
        .errorhandling = "pass") %dopar% {

  fname <- basename(gwas_file)
  pheno_name <- extract_pheno(fname)
  pheno_tag  <- safe_name(pheno_name)

  gw_n <- gwas_sample_df[file == fname, N]
  if (length(gw_n) == 0) return(NULL)

  out_pheno_dir <- file.path(out_root, pheno_tag)
  dir.create(out_pheno_dir, showWarnings = FALSE, recursive = TRUE)

  message(paste(Sys.time(), "[GWAS] Processing:", pheno_name, "N =", gw_n))

  gwas <- tryCatch(fread(gwas_file),
                   error = function(e) {
                     message("GWAS read failed: ", fname, " - ", e$message)
                     NULL
                   })
  if (is.null(gwas)) return(NULL)

  needed_gwas <- c("Chr", "SNP", "bp", "A1", "A2", "Freq", "b", "se", "p")
  if (!all(needed_gwas %in% names(gwas))) {
    message("GWAS missing columns: ",
            paste(setdiff(needed_gwas, names(gwas)), collapse = ", "))
    return(NULL)
  }

  gwas[, maf := ifelse(Freq <= 0.5, Freq, 1 - Freq)]
  setDT(gwas)
  setkey(gwas, SNP)

  for (eqtl_info in eqtl_info_list) {
    period    <- eqtl_info$Period
    tissue    <- eqtl_info$Tissue
    eqtl_file <- eqtl_info$file

    eqtl_n <- eqtl_sample[Period == period & Tissue == tissue, SampleSize]
    if (length(eqtl_n) == 0) next

    message(paste("  [eQTL]", basename(eqtl_file), "N =", eqtl_n))

    eqtl <- tryCatch(fread(eqtl_file),
                     error = function(e) {
                       message("eQTL read failed: ",
                               basename(eqtl_file), " - ", e$message)
                       NULL
                     })
    if (is.null(eqtl)) next

    needed_eqtl <- c("pheno_id", "variant_id", "af",
                     "beta_g1", "beta_se_g1", "pval_g1")
    if (!all(needed_eqtl %in% names(eqtl))) {
      message("eQTL missing columns: ",
              paste(setdiff(needed_eqtl, names(eqtl)), collapse = ", "))
      next
    }

    eqtl <- eqtl[, .(phenotype_id = pheno_id, variant_id, af,
                     slope = beta_g1,
                     slope_se = beta_se_g1,
                     pval_nominal = pval_g1)]
    eqtl[, af := ifelse(af <= 0.5, af, 1 - af)]
    setDT(eqtl)
    setkey(eqtl, variant_id)

    genes <- unique(eqtl$phenotype_id)
    gene_idx <- 0

    for (gene in genes) {
      qtl_sub <- eqtl[phenotype_id == gene]
      if (nrow(qtl_sub) == 0) next

      common <- intersect(qtl_sub$variant_id, gwas$SNP)
      if (length(common) < 2) next

      df <- merge(gwas[J(common)], qtl_sub[J(common)],
                  by.x = "SNP", by.y = "variant_id")
      if (nrow(df) == 0) next

      if (anyNA(df[, .(b, se, p, maf,
                       slope, slope_se, af, pval_nominal)])) next

      d1 <- list(beta = df$b, varbeta = df$se^2, snp = df$SNP,
                 type = "quant",
                 N = rep(gw_n, nrow(df)),
                 MAF = df$maf)
      d2 <- list(beta = df$slope, varbeta = df$slope_se^2, snp = df$SNP,
                 type = "quant",
                 N = rep(eqtl_n, nrow(df)),
                 MAF = df$af)

      res <- tryCatch(coloc.abf(d1, d2), error = function(e) NULL)
      if (is.null(res)) next

      lead_snp <- df$SNP[which.min(df$p)]

      gene_tag   <- safe_name(gene)
      tissue_tag <- safe_name(tissue)
      period_tag <- safe_name(period)

      summary_df <- as.data.frame(t(unlist(res$summary)))
      summary_df$lead_SNP  <- lead_snp
      summary_df$gene      <- gene
      summary_df$tissue    <- tissue
      summary_df$period    <- period
      summary_df$phenotype <- pheno_name

      details_df <- as.data.frame(res$results)
      details_df$gene      <- gene
      details_df$tissue    <- tissue
      details_df$period    <- period
      details_df$phenotype <- pheno_name

      sum_file <- file.path(
        out_pheno_dir,
        paste0(pheno_tag, "__", gene_tag, "__", tissue_tag,
               "__", period_tag, "__summary.txt"))
      det_file <- file.path(
        out_pheno_dir,
        paste0(pheno_tag, "__", gene_tag, "__", tissue_tag,
               "__", period_tag, "__details.txt"))

      fwrite(summary_df, sum_file, sep = "\t")
      fwrite(details_df, det_file, sep = "\t")

      gene_idx <- gene_idx + 1
      if (gene_idx %% 10 == 0) {
        message(paste("    genes done:", gene_idx, "/", length(genes)))
      }
    }
    message(paste("  [eQTL] done:", basename(eqtl_file)))
  }
  message(paste("[GWAS] done:", pheno_name))
}

stopCluster(cl)
log_msg("====== Colocalization completed ======")