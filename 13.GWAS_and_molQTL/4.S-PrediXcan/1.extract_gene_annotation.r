#!/usr/bin/env Rscript
#
# extract_gene_annotation.r
# Extract gene-level annotation from a GTF file and write a clean table.
#
# Arguments:
#   gtf_file   Path to the GTF file
#   out_dir    Output directory (created if missing)
#   out_name   Output filename (default: gene_annotation.txt)
#

suppressPackageStartupMessages({
  library(data.table)
})

#  Parse arguments 
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage: Rscript extract_gene_annotation.R <gtf_file> <out_dir> [out_name]\n")
  quit(status = 1)
}

gtf_file <- args[1]
out_dir  <- args[2]
out_name <- if (length(args) >= 3) args[3] else "gene_annotation.txt"

if (!file.exists(gtf_file)) {
  stop("GTF file not found: ", gtf_file)
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_dir, out_name)

#  Read GTF 
gtf <- fread(gtf_file, header = FALSE, sep = "\t", data.table = FALSE)
colnames(gtf) <- paste0("V", seq_len(ncol(gtf)))

# Keep only gene-level rows
gene_rows <- gtf[gtf$V3 == "gene", ]
if (nrow(gene_rows) == 0) {
  stop("No gene-level records found in ", gtf_file)
}

#  Attribute parser 
get_attr_value <- function(attr_string, key = NULL, is_dbxref_geneid = FALSE) {
  if (is_dbxref_geneid) {
    pattern <- 'db_xref "GeneID:([0-9]+)"'
  } else {
    pattern <- paste0(key, ' "([^"]+)"')
  }
  m <- regmatches(attr_string, regexpr(pattern, attr_string))
  if (length(m) == 0) return(NA_character_)
  sub(pattern, "\\1", m)
}

#  Extract fields 
gene_id      <- vapply(gene_rows$V9, get_attr_value,
                       FUN.VALUE = character(1), is_dbxref_geneid = TRUE)
gene_name    <- vapply(gene_rows$V9, get_attr_value,
                       FUN.VALUE = character(1), key = "gene")
gene_biotype <- vapply(gene_rows$V9, get_attr_value,
                       FUN.VALUE = character(1), key = "gene_biotype")

result <- data.frame(
  chr       = gene_rows$V1,
  gene_id   = gene_id,
  gene_name = gene_name,
  start     = gene_rows$V4,
  end       = gene_rows$V5,
  gene_type = gene_biotype,
  stringsAsFactors = FALSE
)

#  Write output 
write.table(result, out_file, row.names = FALSE, quote = FALSE, sep = "\t")
cat("Gene annotation written to:", out_file, "\n")