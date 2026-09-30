#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(data.table)})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {stop("Usage: Rscript polarize_snps.R <input.tsv> <output.tsv> [fix_thr] [poly_thr] [anc_thr]")}
input_file  <- args[1]
output_file <- args[2]
fix_thr  <- ifelse(length(args) >= 3, as.numeric(args[3]), 0.01)
poly_thr <- ifelse(length(args) >= 4, as.numeric(args[4]), 0.05)
anc_thr  <- ifelse(length(args) >= 5, as.numeric(args[5]), 0.50)
d <- fread(input_file, header = TRUE)
required_cols <- c("Broiler_FileALT_AF", "XJH_FileALT_AF", "RJF_FileALT_AF")
missing_cols <- setdiff(required_cols, names(d))
if (length(missing_cols) > 0) {stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))}

maf <- function(p) pmin(p, 1 - p)
freqC_of_fixed <- function(pop_altfreq, altfreq_C) {ifelse(pop_altfreq > 0.5, altfreq_C, 1 - altfreq_C)}
polarize <- function(fC) {
  out <- rep("ancestral_unsure", length(fC))
  out[fC >= anc_thr]       <- "ancestral_fixed"
  out[fC <= (1 - anc_thr)] <- "derived_fixed"
  out
}

d$maf_A <- maf(d$Broiler_FileALT_AF)
d$maf_B <- maf(d$XJH_FileALT_AF)

d$group <- "other"
d$group[d$maf_A > poly_thr & d$maf_B < fix_thr]  <- "Apoly_Bfix"
d$group[d$maf_A < fix_thr  & d$maf_B > poly_thr] <- "Afix_Bpoly"
d$group[d$maf_A > poly_thr & d$maf_B > poly_thr] <- "Apoly_Bpoly"

fC_A <- freqC_of_fixed(d$Broiler_FileALT_AF, d$RJF_FileALT_AF)
fC_B <- freqC_of_fixed(d$XJH_FileALT_AF, d$RJF_FileALT_AF)

d$freqC_of_fixed <- NA_real_
d$polarized      <- "not_applicable"

isAfix <- d$group == "Afix_Bpoly"
d$freqC_of_fixed[isAfix] <- round(fC_A[isAfix], 4)
d$polarized[isAfix]      <- paste0("A_", polarize(fC_A[isAfix]))

isBfix <- d$group == "Apoly_Bfix"
d$freqC_of_fixed[isBfix] <- round(fC_B[isBfix], 4)
d$polarized[isBfix]      <- paste0("B_", polarize(fC_B[isBfix]))

d$polar_conf <- ifelse(abs(d$freqC_of_fixed - 0.5) >= 0.3, "high", "low")

fwrite(d, output_file, sep = "\t", quote = FALSE)