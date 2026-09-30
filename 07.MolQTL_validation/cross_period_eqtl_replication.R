#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(qvalue)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(
    "Usage: Rscript cross_period_eqtl_replication.R ",
    "<discovery_dir> <validation_dir> <output_file>"
  )
}

discovery_dir <- args[1]
validation_dir <- args[2]
output_file <- args[3]

periods <- c("A", "B", "C", "ALL")

tissues <- c(
  "AF", "CT", "FZ", "GZ", "HC", "JG", "KC", "MC",
  "PZ", "SWM", "SZ", "XJ", "XQN", "XW", "XZ", "ZC"
)

chromosomes <- 1:39


# --------------------------------------------------
# Calculate pi1
# --------------------------------------------------

calc_pi1 <- function(pvalues) {

  pvalues <- pvalues[
    is.finite(pvalues) &
    pvalues >= 0 &
    pvalues <= 1
  ]

  if (length(pvalues) < 2)
    return(NA_real_)

  qobj <- try(
    qvalue(p = pvalues),
    silent = TRUE
  )

  if (inherits(qobj, "try-error")) {

    qobj <- try(
      qvalue(p = c(pvalues, 1)),
      silent = TRUE
    )

    if (inherits(qobj, "try-error"))
      return(NA_real_)
  }

  1 - qobj$pi0
}


# --------------------------------------------------
# Cross-period replication
# --------------------------------------------------

results <- list()
k <- 1

for (discovery_period in periods) {

  for (validation_period in setdiff(periods, discovery_period)) {

    for (tissue in tissues) {

      cat(
        "Discovery:", discovery_period,
        "| Validation:", validation_period,
        "| Tissue:", tissue, "\n"
      )

      discovery_file <- file.path(
        discovery_dir,
        discovery_period,
        paste0(
          discovery_period,
          "_",
          tissue,
          "_LMM.cis_qtl.txt"
        )
      )

      if (!file.exists(discovery_file))
        next

      discovery <- fread(discovery_file)

      required_cols <- c(
        "pheno_id",
        "variant_id",
        "beta_g1",
        "beta_se_g1",
        "pval_g1"
      )

      if (!all(required_cols %in% colnames(discovery)))
        next

      validation_p <- numeric()
      discovery_z <- numeric()
      validation_z <- numeric()

      for (chr in chromosomes) {

        validation_file <- file.path(
          validation_dir,
          validation_period,
          tissue,
          paste0(
            validation_period,
            "_",
            tissue,
            "_LMM.cis_qtl_pairs.",
            chr,
            ".txt"
          )
        )

        if (!file.exists(validation_file))
          next

        validation <- fread(validation_file)

        if (!all(required_cols %in% colnames(validation)))
          next

        matched <- merge(
          discovery,
          validation,
          by = c("pheno_id", "variant_id"),
          suffixes = c(".discovery", ".validation")
        )

        if (nrow(matched) == 0)
          next

        validation_p <- c(
          validation_p,
          matched$pval_g1.validation
        )

        discovery_z <- c(
          discovery_z,
          matched$beta_g1.discovery /
            matched$beta_se_g1.discovery
        )

        validation_z <- c(
          validation_z,
          matched$beta_g1.validation /
            matched$beta_se_g1.validation
        )
      }

      if (length(validation_p) == 0)
        next

      keep <- (
        is.finite(discovery_z) &
        is.finite(validation_z)
      )

      pi1 <- calc_pi1(validation_p)

      rho <- if (sum(keep) >= 2) {
        cor(
          discovery_z[keep],
          validation_z[keep],
          method = "spearman"
        )
      } else {
        NA_real_
      }

      results[[k]] <- data.table(
        DiscoveryPeriod = discovery_period,
        ValidationPeriod = validation_period,
        Tissue = tissue,
        MatchedPairs = length(validation_p),
        Pi1 = pi1,
        SpearmanCorrelation = rho
      )

      k <- k + 1
    }
  }
}


# --------------------------------------------------
# Output
# --------------------------------------------------

if (length(results) == 0)
  stop("No matched eQTL pairs found.")

results <- rbindlist(results)

fwrite(
  results,
  output_file,
  sep = "\t"
)

cat("Finished.\n")
cat("Comparisons:", nrow(results), "\n")
cat("Output:", output_file, "\n")