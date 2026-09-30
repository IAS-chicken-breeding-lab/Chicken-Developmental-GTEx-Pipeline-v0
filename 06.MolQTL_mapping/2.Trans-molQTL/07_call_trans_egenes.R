#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: 07_call_trans_egenes.R INPUT OUT_PREFIX FDR EFFECTIVE_TRANS_TESTS")
}

input_file <- args[1]
out_prefix <- args[2]
fdr_threshold <- as.numeric(args[3])
effective_tests <- as.numeric(args[4])

if (!file.exists(input_file)) stop("Missing input: ", input_file)
if (!is.finite(fdr_threshold) || fdr_threshold <= 0 || fdr_threshold >= 1) {
  stop("FDR must be between 0 and 1")
}
if (!is.finite(effective_tests) || effective_tests < 1) stop("Invalid EFFECTIVE_TRANS_TESTS")

trans <- fread(input_file)
required <- c("pheno_id", "variant_id", "pval_g1")
missing <- setdiff(required, names(trans))
if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))

trans[, pval_g1 := suppressWarnings(as.numeric(pval_g1))]
if (any(!is.finite(trans$pval_g1)) || any(trans$pval_g1 < 0 | trans$pval_g1 > 1)) {
  stop("pval_g1 contains non-numeric or out-of-range values")
}
if (!"tissue" %in% names(trans)) {
  trans[, tissue := sub("\\..*$", "", basename(input_file))]
}

# Reproduce the correction used in the original project:
# 1) multiply nominal P by the fixed effective number of trans tests;
# 2) BH-correct those values within each phenotype gene;
# 3) take one deterministic lead variant per gene and BH-correct across genes.
trans[, pval_multiple := pmin(pval_g1 * effective_tests, 1)]
trans[, pval_adj_BH_gene := p.adjust(pval_multiple, method = "BH"), by = pheno_id]
setorder(trans, pheno_id, pval_multiple, variant_id)
lead <- trans[, .SD[1], by = pheno_id]
lead[, qval := p.adjust(pval_multiple, method = "BH")]
lead[, is_eGene := qval < fdr_threshold & pval_adj_BH_gene < fdr_threshold]

gene_status <- lead[, .(pheno_id, qval, is_eGene)]
trans_sig <- merge(trans, gene_status, by = "pheno_id", all.x = TRUE, sort = FALSE)
trans_sig <- trans_sig[is_eGene == TRUE & pval_adj_BH_gene < fdr_threshold]

fwrite(lead, paste0(out_prefix, ".txt.gz"), sep = "\t")
fwrite(trans_sig, paste0(out_prefix, ".sig.txt.gz"), sep = "\t")

cat("Genes tested:", nrow(lead), "\n")
cat("trans-eGenes:", sum(lead$is_eGene), "\n")
cat("Significant pairs:", nrow(trans_sig), "\n")

