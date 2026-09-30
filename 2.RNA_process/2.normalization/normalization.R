#!/usr/bin/env Rscript

# ==============================================================================
# featureCounts_to_TPM.R
#
# Description:
#   Merge single-sample featureCounts files, convert raw counts to TPM,
#   perform basic gene/sample expression filtering, and export matched
#   count and TPM matrices.
#
# Input:
#   A directory containing single-sample featureCounts output files.
#
# Default filtering:
#   1. Remove genes with TPM < 0.1 in all samples.
#   2. Remove samples in which fewer than 20% of retained genes have TPM >= 0.1.
#
# Output:
#   <prefix>.counts.tsv
#   <prefix>.tpm.tsv
#   <prefix>.filtered_genes.txt
#   <prefix>.filtered_samples.txt
#   <prefix>.qc_summary.tsv
#
# Example:
#   Rscript featureCounts_to_TPM.R \
#       --input_path /path/to/featureCounts \
#       --pattern '[.]tsv$' \
#       --output_path /path/to/output \
#       --prefix expression
#
# Dependencies:
#   R >= 4.0
#   data.table
#   argparser
#
# ==============================================================================

options(warn = -1)

suppressPackageStartupMessages({
    library(data.table)
    library(argparser)
})


# ==============================================================================
# Arguments
# ==============================================================================

p <- arg_parser(
    "Merge featureCounts files and calculate TPM"
)

p <- add_argument(
    p,
    "--input_path",
    help = "Directory containing featureCounts files",
    type = "character",
    default = "./"
)

p <- add_argument(
    p,
    "--pattern",
    help = "Regular expression used to select input files",
    type = "character",
    default = "\\.tsv$"
)

p <- add_argument(
    p,
    "--output_path",
    help = "Output directory",
    type = "character",
    default = "./"
)

p <- add_argument(
    p,
    "--prefix",
    help = "Prefix for output files",
    type = "character",
    default = "expression"
)

p <- add_argument(
    p,
    "--min_tpm",
    help = "TPM threshold defining an expressed gene",
    type = "numeric",
    default = 0.1
)

p <- add_argument(
    p,
    "--min_expressed_fraction",
    help = "Minimum fraction of expressed genes required for a sample",
    type = "numeric",
    default = 0.2
)

argv <- parse_args(p)

input_path <- argv$input_path
pattern <- argv$pattern
output_path <- argv$output_path
prefix <- argv$prefix
min_tpm <- argv$min_tpm
min_expressed_fraction <- argv$min_expressed_fraction


# ==============================================================================
# Basic checks
# ==============================================================================

if (!dir.exists(input_path)) {
    stop("Input directory does not exist: ", input_path)
}

if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
}

if (min_tpm < 0) {
    stop("--min_tpm must be >= 0")
}

if (min_expressed_fraction <= 0 ||
    min_expressed_fraction > 1) {
    stop("--min_expressed_fraction must be between 0 and 1")
}


# ==============================================================================
# Find input files
# ==============================================================================

files <- list.files(
    path = input_path,
    pattern = pattern,
    full.names = TRUE
)

files <- sort(files)

if (length(files) == 0) {
    stop(
        "No files found in ",
        input_path,
        " using pattern: ",
        pattern
    )
}

cat("============================================================\n")
cat("featureCounts to TPM\n")
cat("============================================================\n")
cat("Input directory:      ", input_path, "\n")
cat("Number of files:      ", length(files), "\n")
cat("TPM threshold:        ", min_tpm, "\n")
cat("Sample gene fraction: ", min_expressed_fraction, "\n")
cat("Output directory:     ", output_path, "\n")
cat("============================================================\n\n")


# ==============================================================================
# Read one featureCounts file
# ==============================================================================

read_featurecounts <- function(file) {

    # featureCounts output often contains a comment line before the header.
    # skip = "Geneid" starts reading from the actual table header.
    x <- fread(
        file,
        skip = "Geneid",
        data.table = FALSE,
        check.names = FALSE
    )

    required_columns <- c(
        "Geneid",
        "Chr",
        "Start",
        "End",
        "Strand",
        "Length"
    )

    if (!all(required_columns %in% colnames(x))) {
        stop(
            "Invalid featureCounts format: ",
            basename(file)
        )
    }

    count_columns <- setdiff(
        colnames(x),
        required_columns
    )

    if (length(count_columns) != 1) {
        stop(
            basename(file),
            " contains ",
            length(count_columns),
            " count columns. ",
            "This script expects one sample per file."
        )
    }

    if (anyDuplicated(x$Geneid)) {
        stop(
            "Duplicated Geneid detected in ",
            basename(file)
        )
    }

    data.frame(
        Geneid = as.character(x$Geneid),
        Length = as.numeric(x$Length),
        Count = as.numeric(x[[count_columns]]),
        stringsAsFactors = FALSE
    )
}


# ==============================================================================
# Read reference file
# ==============================================================================

cat("Reading featureCounts files...\n")

ref <- read_featurecounts(files[1])

genes <- ref$Geneid
gene_length <- ref$Length

if (any(is.na(gene_length)) || any(gene_length <= 0)) {
    stop("Invalid gene lengths detected.")
}

sample_names <- sub(
    "\\.[^.]+$",
    "",
    basename(files)
)

if (anyDuplicated(sample_names)) {
    stop(
        "Duplicated sample names detected after removing file extensions."
    )
}


# ==============================================================================
# Allocate count matrix
# ==============================================================================

count_mat <- matrix(
    NA_real_,
    nrow = length(genes),
    ncol = length(files),
    dimnames = list(
        genes,
        sample_names
    )
)

count_mat[, 1] <- ref$Count

cat(
    sprintf(
        "[%d/%d] %s\n",
        1,
        length(files),
        basename(files[1])
    )
)


# ==============================================================================
# Read remaining files
# ==============================================================================

if (length(files) > 1) {

    for (i in 2:length(files)) {

        x <- read_featurecounts(files[i])

        # All files should contain the same genes
        if (!setequal(genes, x$Geneid)) {
            stop(
                "Gene sets differ between files:\n",
                basename(files[1]),
                "\n",
                basename(files[i])
            )
        }

        # Match order to reference
        idx <- match(genes, x$Geneid)

        # Check gene lengths
        if (!all(gene_length == x$Length[idx])) {
            stop(
                "Gene lengths differ in file: ",
                basename(files[i])
            )
        }

        count_mat[, i] <- x$Count[idx]

        cat(
            sprintf(
                "[%d/%d] %s\n",
                i,
                length(files),
                basename(files[i])
            )
        )
    }
}


# ==============================================================================
# Count matrix checks
# ==============================================================================

if (anyNA(count_mat)) {
    stop("NA values detected in count matrix.")
}

if (any(count_mat < 0)) {
    stop("Negative counts detected.")
}

cat("\nCount matrix successfully generated.\n")
cat("Genes:   ", nrow(count_mat), "\n")
cat("Samples: ", ncol(count_mat), "\n\n")


# ==============================================================================
# Convert counts to TPM
#
# TPM_i =
#
#      count_i / length_i
# ----------------------------- x 1,000,000
# sum(count_j / length_j)
#
# Gene length is supplied by featureCounts in base pairs.
# ==============================================================================

counts_to_tpm <- function(counts, gene_length) {

    gene_length_kb <- gene_length / 1000

    # Reads/fragments per kilobase
    rpk <- sweep(
        counts,
        1,
        gene_length_kb,
        FUN = "/"
    )

    scaling_factor <- colSums(rpk)

    if (any(scaling_factor <= 0)) {
        stop(
            "At least one sample has zero total RPK."
        )
    }

    tpm <- sweep(
        rpk,
        2,
        scaling_factor,
        FUN = "/"
    ) * 1e6

    return(tpm)
}


cat("Calculating TPM...\n")

tpm <- counts_to_tpm(
    count_mat,
    gene_length
)


# ==============================================================================
# TPM sanity check
# ==============================================================================

tpm_sum <- colSums(tpm)

if (any(abs(tpm_sum - 1e6) > 1)) {
    warning(
        "Some samples do not sum to approximately 1,000,000 TPM."
    )
}

cat(
    "TPM calculation completed.\n",
    "TPM column sums before filtering: ",
    sprintf(
        "%.2f - %.2f",
        min(tpm_sum),
        max(tpm_sum)
    ),
    "\n\n",
    sep = ""
)


# ==============================================================================
# Gene filtering
#
# Remove genes with TPM < min_tpm in every sample.
#
# Equivalent to:
#   keep gene if TPM >= min_tpm in at least one sample.
# ==============================================================================

gene_keep <- rowSums(
    tpm >= min_tpm
) > 0

filtered_genes <- rownames(tpm)[
    !gene_keep
]

count_mat <- count_mat[
    gene_keep,
    ,
    drop = FALSE
]

tpm <- tpm[
    gene_keep,
    ,
    drop = FALSE
]

cat(
    length(filtered_genes),
    " genes were removed.\n"
)

cat(
    nrow(tpm),
    " genes remained.\n\n"
)


# ==============================================================================
# Sample filtering
#
# Remove samples where fewer than min_expressed_fraction of retained genes
# have TPM >= min_tpm.
# ==============================================================================

expressed_gene_fraction <- colMeans(
    tpm >= min_tpm
)

sample_keep <- expressed_gene_fraction >=
    min_expressed_fraction

filtered_samples <- colnames(tpm)[
    !sample_keep
]

count_mat <- count_mat[
    ,
    sample_keep,
    drop = FALSE
]

tpm <- tpm[
    ,
    sample_keep,
    drop = FALSE
]

cat(
    length(filtered_samples),
    " samples were removed.\n"
)

cat(
    ncol(tpm),
    " samples remained.\n\n"
)


# ==============================================================================
# Output
# ==============================================================================

count_file <- file.path(
    output_path,
    paste0(prefix, ".counts.tsv")
)

tpm_file <- file.path(
    output_path,
    paste0(prefix, ".tpm.tsv")
)

filtered_gene_file <- file.path(
    output_path,
    paste0(prefix, ".filtered_genes.txt")
)

filtered_sample_file <- file.path(
    output_path,
    paste0(prefix, ".filtered_samples.txt")
)

qc_file <- file.path(
    output_path,
    paste0(prefix, ".qc_summary.tsv")
)


# Counts
count_out <- as.data.table(
    count_mat,
    keep.rownames = "gene_id"
)

fwrite(
    count_out,
    count_file,
    sep = "\t",
    quote = FALSE
)


# TPM
tpm_out <- as.data.table(
    tpm,
    keep.rownames = "gene_id"
)

fwrite(
    tpm_out,
    tpm_file,
    sep = "\t",
    quote = FALSE
)


# Filtered genes
writeLines(
    filtered_genes,
    filtered_gene_file
)


# Filtered samples
writeLines(
    filtered_samples,
    filtered_sample_file
)


# ==============================================================================
# QC summary
# ==============================================================================

qc_summary <- data.frame(
    Metric = c(
        "Input_files",
        "Input_genes",
        "Genes_after_filtering",
        "Filtered_genes",
        "Input_samples",
        "Samples_after_filtering",
        "Filtered_samples",
        "TPM_threshold",
        "Minimum_expressed_gene_fraction"
    ),

    Value = c(
        length(files),
        length(genes),
        nrow(tpm),
        length(filtered_genes),
        length(files),
        ncol(tpm),
        length(filtered_samples),
        min_tpm,
        min_expressed_fraction
    )
)

fwrite(
    qc_summary,
    qc_file,
    sep = "\t",
    quote = FALSE
)


# ==============================================================================
# Finish
# ==============================================================================

cat("============================================================\n")
cat("Completed successfully\n")
cat("============================================================\n")
cat("Counts:           ", count_file, "\n")
cat("TPM:              ", tpm_file, "\n")
cat("Filtered genes:   ", filtered_gene_file, "\n")
cat("Filtered samples: ", filtered_sample_file, "\n")
cat("QC summary:       ", qc_file, "\n")
cat("============================================================\n")