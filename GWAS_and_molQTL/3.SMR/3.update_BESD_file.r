#!/usr/bin/env Rscript
#
# Update SMR BESD .esi and .epi files with correct SNP frequency and gene annotation.
#
# Usage:
#   Rscript update_besd.R <tss_file> <gtf_file> <frq_dir> <besd_dir> <output_dir> <tissue>
#
# Arguments:
#   tss_file    TSS annotation file (chr, start, strand, gene_id, gene_name)
#   gtf_file    GTF file (currently unused, kept for compatibility)
#   frq_dir     Directory containing <tissue>_freq.frq
#   besd_dir    Directory containing original <tissue>_eGene_ciseQTL.esi and .epi
#   output_dir  Directory for updated .esi and .epi files
#   tissue      Tissue identifier (e.g., A_AF)
#


suppressMessages(library(data.table))

`%&%` <- function(a, b) paste0(a, b)

#  Parse arguments 
ARGS <- commandArgs(trailingOnly = TRUE)

if (length(ARGS) < 6) {
  cat("Usage: Rscript update_besd.R <tss_file> <gtf_file> <frq_dir> <besd_dir> <output_dir> <tissue>\n")
  quit(status = 1)
}

file_tss      <- ARGS[1]
file_gtf      <- ARGS[2]
dir_frq       <- ARGS[3]
dir_BESD_file <- ARGS[4]
output_dir    <- ARGS[5]
tis           <- ARGS[6]

message("Processing tissue: ", tis)

#  Read annotation 
if (!file.exists(file_tss)) stop("Annotation file not found: ", file_tss)
annot <- fread(file_tss)
required_annot <- c("chr", "start", "strand", "gene_id", "gene_name")
if (!all(required_annot %in% names(annot))) {
  stop("Annotation file must contain columns: ", paste(required_annot, collapse = ", "))
}

#  Read frequency file 
frq_file <- file.path(dir_frq, paste0(tis, "_freq.frq"))
if (!file.exists(frq_file)) stop("Frequency file not found: ", frq_file)
frq <- fread(frq_file)
required_frq <- c("CHR", "SNP", "A1", "A2", "MAF")
if (!all(required_frq %in% names(frq))) {
  stop("Frequency file must contain columns: ", paste(required_frq, collapse = ", "))
}

#  Read old ESI / EPI 
esi_file <- file.path(dir_BESD_file, paste0(tis, "_eGene_ciseQTL.esi"))
epi_file <- file.path(dir_BESD_file, paste0(tis, "_eGene_ciseQTL.epi"))

if (!file.exists(esi_file)) stop("ESI file not found: ", esi_file)
if (!file.exists(epi_file)) stop("EPI file not found: ", epi_file)

esi <- fread(esi_file)
epi <- fread(epi_file)

#  Update ESI 
message("Updating ESI ...")

matIdx <- match(esi$V2, frq$SNP)
if (any(is.na(matIdx))) {
  warning(sum(is.na(matIdx)), " SNPs in ESI not found in frequency file.")
}

esi$V1 <- frq$CHR[matIdx]
esi$V5 <- frq$A1[matIdx]
esi$V6 <- frq$A2[matIdx]
esi$V7 <- frq$MAF[matIdx]
esi$V4 <- gsub("^[0-9]+_|_(A|T|C|G)$", "", esi$V2)

output_esi_file <- file.path(output_dir, paste0(tis, "_eGene_ciseQTL.esi"))
fwrite(esi, output_esi_file, col.names = FALSE, sep = "\t")
message("ESI updated: ", output_esi_file)

#  Update EPI 
message("Updating EPI ...")

annotIdx <- match(epi$V2, annot$gene_name)
if (all(is.na(annotIdx))) {
  stop("None of the gene names in EPI matched annotation file.")
} else if (any(is.na(annotIdx))) {
  warning(sum(is.na(annotIdx)), " genes in EPI not found in annotation.")
}

epi$V1 <- annot$chr[annotIdx]
epi$V2 <- annot$gene_id[annotIdx]
epi$V4 <- annot$start[annotIdx]
epi$V5 <- annot$gene_name[annotIdx]
epi$V6 <- annot$strand[annotIdx]

output_epi_file <- file.path(output_dir, paste0(tis, "_eGene_ciseQTL.epi"))
fwrite(epi, output_epi_file, col.names = FALSE, sep = "\t")
message("EPI updated: ", output_epi_file)

message("Finished processing: ", tis)