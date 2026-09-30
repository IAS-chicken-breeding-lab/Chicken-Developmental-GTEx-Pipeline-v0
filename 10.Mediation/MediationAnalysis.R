#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(mediation)
  library(data.table)
  library(argparse)
})

# --------------------------------------------------
# Arguments
# --------------------------------------------------

parser <- ArgumentParser()

parser$add_argument("--period", required = TRUE)
parser$add_argument("--tissue", required = TRUE)
parser$add_argument("--gene", required = TRUE)

parser$add_argument("--cell")
parser$add_argument("--eigengene")
parser$add_argument("--tf")
parser$add_argument("--hormone")

parser$add_argument("--snp", required = TRUE)
parser$add_argument("--genotype_file", required = TRUE)
parser$add_argument("--cov_file", required = TRUE)
parser$add_argument("--output_dir", required = TRUE)

args <- parser$parse_args()

mediators <- c(
  cell = args$cell,
  eigengene = args$eigengene,
  tf = args$tf,
  hormone = args$hormone
)

provided <- which(!vapply(
  mediators,
  is.null,
  logical(1)
))

if (length(provided) != 1) {
  stop("Exactly one mediator must be specified.")
}

mediator_type <- names(provided)
mediator_name <- mediators[[provided]]

# --------------------------------------------------
# Helper
# --------------------------------------------------

INT <- function(x) {

  x <- as.numeric(as.character(x))

  out <- rep(NA_real_, length(x))

  idx <- which(!is.na(x))

  r <- rank(
    x[idx],
    ties.method = "average"
  )

  out[idx] <- qnorm(
    (r - 0.5) / length(idx)
  )

  out
}

# --------------------------------------------------
# Read data
# --------------------------------------------------

dat <- fread(
  args$genotype_file,
  data.table = FALSE,
  check.names = FALSE
)

if (!(args$gene %in% colnames(dat))) {
  stop("Gene not found: ", args$gene)
}

if (!(mediator_name %in% colnames(dat))) {
  stop("Mediator not found: ", mediator_name)
}

# SNP
if (args$snp %in% colnames(dat)) {

  snp_col <- args$snp

} else {

  snp_col <- grep(
    args$snp,
    colnames(dat),
    value = TRUE,
    fixed = TRUE
  )

  if (length(snp_col) != 1) {
    stop("SNP cannot be uniquely identified: ", args$snp)
  }
}

# --------------------------------------------------
# Covariates
# --------------------------------------------------

cov <- fread(
  args$cov_file,
  data.table = FALSE,
  check.names = FALSE
)

cov <- as.data.frame(t(cov))

colnames(cov) <- cov[1, ]
cov <- cov[-1, , drop = FALSE]

cov$ID <- rownames(cov)

if (!all(dat$ID %in% cov$ID)) {
  stop("Some genotype IDs are absent from the covariate file.")
}

dat <- merge(
  dat,
  cov,
  by = "ID"
)

# Do not include raw Sex twice
covariate_columns <- setdiff(
  colnames(cov),
  c("ID", "Sex")
)

# --------------------------------------------------
# Variables
# --------------------------------------------------

numeric_columns <- unique(c(
  args$gene,
  snp_col,
  mediator_name,
  "Sex",
  covariate_columns
))

for (x in numeric_columns) {

  if (x %in% colnames(dat)) {
    dat[[x]] <- suppressWarnings(
      as.numeric(as.character(dat[[x]]))
    )
  }
}

dat$Y <- dat[[args$gene]]

# Mediator
dat$mediator_value <- dat[[mediator_name]]

# Cell abundance QC
if (mediator_type == "cell") {

  if (all(dat$mediator_value == 0, na.rm = TRUE)) {
    stop("Cell abundance is zero in all samples.")
  }

  mean1 <- mean(
    dat$mediator_value[dat$Sex == 1],
    na.rm = TRUE
  )

  mean2 <- mean(
    dat$mediator_value[dat$Sex == 2],
    na.rm = TRUE
  )

  if (mean1 <= 0.025 && mean2 <= 0.025) {
    stop("Cell abundance <= 0.025 in both sex groups.")
  }
}

# INT for hormone and cell mediator
if (mediator_type %in% c("hormone", "cell")) {

  dat$mediator_value <- INT(
    dat$mediator_value
  )

  dat$mediator_value <- as.numeric(
    scale(
      dat$mediator_value,
      center = TRUE,
      scale = FALSE
    )
  )
}

# SNP
dat$G <- INT(
  dat[[snp_col]]
)

dat$G <- as.numeric(
  scale(
    dat$G,
    center = TRUE,
    scale = FALSE
  )
)

# Sex
dat$sex_centered <- INT(
  dat$Sex
)

dat$sex_centered <- as.numeric(
  scale(
    dat$sex_centered,
    center = TRUE,
    scale = FALSE
  )
)

# Interaction terms
dat$Sex_G <- (
  dat$sex_centered *
  dat$G
)

dat$G_mediator <- (
  dat$G *
  dat$mediator_value
)

# --------------------------------------------------
# Complete cases
# --------------------------------------------------

model_vars <- unique(c(
  "Y",
  "sex_centered",
  "G",
  "Sex_G",
  "mediator_value",
  "G_mediator",
  covariate_columns
))

model_vars <- model_vars[
  model_vars %in% colnames(dat)
]

dat <- dat[
  complete.cases(dat[, model_vars, drop = FALSE]),
  ,
  drop = FALSE
]

# --------------------------------------------------
# Models
# --------------------------------------------------

covar_string <- ""

if (length(covariate_columns) > 0) {

  covar_string <- paste0(
    " + ",
    paste(
      covariate_columns,
      collapse = " + "
    )
  )
}

mediator_formula <- as.formula(
  paste0(
    "G_mediator ~ ",
    "sex_centered + G + Sex_G + mediator_value",
    covar_string
  )
)

outcome_formula <- as.formula(
  paste0(
    "Y ~ ",
    "sex_centered + G + Sex_G + ",
    "mediator_value + G_mediator",
    covar_string
  )
)

mediator_model <- lm(
  mediator_formula,
  data = dat
)

outcome_model <- lm(
  outcome_formula,
  data = dat
)

# --------------------------------------------------
# Mediation
# --------------------------------------------------

result <- mediate(
  mediator_model,
  outcome_model,
  treat = "Sex_G",
  mediator = "G_mediator",
  boot = TRUE,
  sims = 999
)

# --------------------------------------------------
# Output
# --------------------------------------------------

dir.create(
  args$output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

output_file <- file.path(
  args$output_dir,
  paste0(
    args$period, "_",
    args$tissue,
    "_mediation_results_",
    args$gene, "_",
    mediator_name, "_",
    args$snp,
    ".txt"
  )
)

sink(output_file)

cat("Mediation Analysis Results\n")
cat("===========================\n")
cat("Period:", args$period, "\n")
cat("Tissue:", args$tissue, "\n")
cat("Gene:", args$gene, "\n")
cat("Mediator type:", mediator_type, "\n")
cat("Mediator:", mediator_name, "\n")
cat("SNP:", args$snp, "\n")
cat("N:", nrow(dat), "\n\n")

print(summary(result))

sink()

cat("Output:", output_file, "\n")