#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(WGCNA)
  library(data.table)
})

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(
    "Usage: Rscript WGCNA_module_trait_analysis.R ",
    "<expression.tsv> <SamplesInfo.txt> <output_dir> ",
    "[threads] [top_var_genes] [maxBlockSize]"
  )
}

expr_file <- args[1]
meta_file <- args[2]
outdir <- args[3]

nThreads <- ifelse(length(args) >= 4, as.integer(args[4]), 32)
top_var_genes <- ifelse(length(args) >= 5, as.integer(args[5]), 18000)
maxBlockSize_use <- ifelse(length(args) >= 6, as.integer(args[6]), 6000)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

Sys.setenv(
  OPENBLAS_NUM_THREADS = 1,
  OMP_NUM_THREADS = 1,
  MKL_NUM_THREADS = 1
)

enableWGCNAThreads(nThreads = nThreads)

# --------------------------------------------------
# Read expression matrix
# --------------------------------------------------

expr <- fread(
  expr_file,
  data.table = FALSE,
  check.names = FALSE
)

gene_ids <- expr[[1]]

expr <- expr[, -1, drop = FALSE]
rownames(expr) <- gene_ids

expr <- as.matrix(expr)
storage.mode(expr) <- "numeric"

# Remove duplicated genes
expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]

# log2(TPM + 1)
expr <- log2(expr + 1)

# --------------------------------------------------
# Select highly variable genes
# --------------------------------------------------

gene_var <- apply(expr, 1, var, na.rm = TRUE)
gene_var[is.na(gene_var)] <- 0

var_table <- data.frame(
  Gene = rownames(expr),
  Variance = gene_var
)

var_table <- var_table[
  order(var_table$Variance, decreasing = TRUE),
]

write.table(
  var_table,
  file.path(outdir, "GeneVariance_AllGenes.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

keep_n <- min(top_var_genes, nrow(var_table))
keep_genes <- var_table$Gene[seq_len(keep_n)]

expr <- expr[keep_genes, , drop = FALSE]

write.table(
  data.frame(Gene = keep_genes),
  file.path(
    outdir,
    paste0("GeneVariance_Top", keep_n, "Genes.txt")
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# WGCNA: samples × genes
datExpr <- as.data.frame(t(expr))

rm(expr)
gc()

# --------------------------------------------------
# Read metadata
# --------------------------------------------------

meta <- fread(
  meta_file,
  data.table = FALSE,
  check.names = FALSE
)

required_cols <- c(
  "BioSample",
  "Tissue",
  "Sex",
  "Period"
)

missing_cols <- setdiff(required_cols, colnames(meta))

if (length(missing_cols) > 0) {
  stop(
    "Missing columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# Match samples
common_samples <- intersect(
  rownames(datExpr),
  meta$BioSample
)

if (length(common_samples) == 0) {
  stop("No matched samples.")
}

datExpr <- datExpr[
  common_samples,
  ,
  drop = FALSE
]

meta <- meta[
  match(common_samples, meta$BioSample),
  ,
  drop = FALSE
]

# --------------------------------------------------
# Check samples and genes
# --------------------------------------------------

gsg <- goodSamplesGenes(
  datExpr,
  verbose = 3
)

if (!gsg$allOK) {

  datExpr <- datExpr[
    gsg$goodSamples,
    gsg$goodGenes,
    drop = FALSE
  ]

  meta <- meta[
    gsg$goodSamples,
    ,
    drop = FALSE
  ]
}

# --------------------------------------------------
# Tissue × Period × Sex traits
# --------------------------------------------------

meta$Group <- paste(
  meta$Tissue,
  meta$Period,
  meta$Sex,
  sep = "_"
)

traitData <- model.matrix(
  ~ 0 + Group,
  data = meta
)

colnames(traitData) <- sub(
  "^Group",
  "",
  colnames(traitData)
)

rownames(traitData) <- meta$BioSample

write.table(
  traitData,
  file.path(outdir, "Trait_Groups.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

write.table(
  meta[, c(
    "BioSample",
    "Tissue",
    "Sex",
    "Period",
    "Group"
  )],
  file.path(outdir, "Sample_Group_Info.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# --------------------------------------------------
# Soft threshold
# --------------------------------------------------

powers <- 1:20

sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  networkType = "signed",
  corFnc = "bicor",
  verbose = 5
)

fit <- sft$fitIndices

write.table(
  fit,
  file.path(outdir, "SoftThreshold.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

pdf(
  file.path(outdir, "SoftThreshold.pdf"),
  width = 12,
  height = 6
)

par(mfrow = c(1, 2))

plot(
  fit[, 1],
  -sign(fit[, 3]) * fit[, 2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit (R²)",
  type = "n"
)

text(
  fit[, 1],
  -sign(fit[, 3]) * fit[, 2],
  labels = fit[, 1]
)

abline(
  h = 0.8,
  lty = 2
)

plot(
  fit[, 1],
  fit[, 5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  type = "n"
)

text(
  fit[, 1],
  fit[, 5],
  labels = fit[, 1]
)

dev.off()

# Select power
candidate <- fit[fit[, 2] >= 0.8, , drop = FALSE]

if (nrow(candidate) > 0) {

  softPower <- candidate[1, 1]

} else {

  candidate <- fit[
    fit[, 2] >= 0.75,
    ,
    drop = FALSE
  ]

  if (nrow(candidate) > 0) {
    softPower <- candidate[1, 1]
  } else {
    softPower <- fit[
      which.max(fit[, 2]),
      1
    ]
  }
}

write.table(
  data.frame(SoftPower = softPower),
  file.path(outdir, "SelectedPower.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# --------------------------------------------------
# WGCNA network
# --------------------------------------------------

net <- blockwiseModules(
  datExpr,
  power = softPower,
  corType = "bicor",
  quickCor = 2,
  networkType = "signed",
  TOMType = "signed",
  maxBlockSize = maxBlockSize_use,
  minModuleSize = 30,
  reassignThreshold = 0,
  deepSplit = 2,
  mergeCutHeight = 0.25,
  numericLabels = FALSE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 5
)

moduleColors <- net$colors

module_gene <- data.frame(
  Gene = colnames(datExpr),
  Module = moduleColors
)

write.table(
  module_gene,
  file.path(outdir, "ModuleGeneMembership.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# --------------------------------------------------
# Module dendrograms
# --------------------------------------------------

pdf(
  file.path(outdir, "ModuleDendrogram.pdf"),
  width = 16,
  height = 8
)

for (i in seq_along(net$dendrograms)) {

  plotDendroAndColors(
    net$dendrograms[[i]],
    moduleColors[
      net$blockGenes[[i]]
    ],
    paste("Module Colors - Block", i),
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05
  )
}

dev.off()

# --------------------------------------------------
# Module size
# --------------------------------------------------

module_size <- sort(
  table(moduleColors),
  decreasing = TRUE
)

write.table(
  module_size,
  file.path(outdir, "ModuleSize.txt"),
  sep = "\t",
  quote = FALSE,
  col.names = FALSE
)

# --------------------------------------------------
# Module eigengenes
# --------------------------------------------------

MEs <- orderMEs(net$MEs)

write.table(
  MEs,
  file.path(outdir, "ModuleEigengenes.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

# --------------------------------------------------
# Module-trait association
# --------------------------------------------------

moduleTraitCor <- cor(
  MEs,
  traitData,
  use = "pairwise.complete.obs"
)

moduleTraitPvalue <- corPvalueStudent(
  moduleTraitCor,
  nrow(datExpr)
)

write.table(
  moduleTraitCor,
  file.path(outdir, "ModuleTraitCor.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

write.table(
  moduleTraitPvalue,
  file.path(outdir, "ModuleTraitPvalue.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

# --------------------------------------------------
# All module-trait pairs
# --------------------------------------------------

cor_df <- as.data.frame(
  as.table(moduleTraitCor)
)

p_df <- as.data.frame(
  as.table(moduleTraitPvalue)
)

colnames(cor_df) <- c(
  "Module",
  "Trait",
  "Correlation"
)

colnames(p_df) <- c(
  "Module",
  "Trait",
  "Pvalue"
)

pair_df <- merge(
  cor_df,
  p_df,
  by = c("Module", "Trait")
)

pair_df$FDR <- p.adjust(
  pair_df$Pvalue,
  method = "BH"
)

pair_df <- pair_df[
  order(
    pair_df$FDR,
    -abs(pair_df$Correlation)
  ),
]

write.table(
  pair_df,
  file.path(outdir, "ModuleTrait_AllPairs.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

pair_sig <- pair_df[
  pair_df$FDR < 0.05 &
    abs(pair_df$Correlation) >= 0.3,
]

write.table(
  pair_sig,
  file.path(
    outdir,
    "ModuleTrait_Significant_FDR005_absCor03.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# --------------------------------------------------
# Module-trait heatmap
# --------------------------------------------------

textMatrix <- paste(
  signif(moduleTraitCor, 2),
  "\n(",
  signif(moduleTraitPvalue, 1),
  ")",
  sep = ""
)

dim(textMatrix) <- dim(moduleTraitCor)

pdf(
  file.path(outdir, "ModuleTraitHeatmap.pdf"),
  width = 55,
  height = 16
)

labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(traitData),
  yLabels = colnames(MEs),
  ySymbols = colnames(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.25,
  zlim = c(-1, 1),
  main = "Module-Trait Relationships"
)

dev.off()

# --------------------------------------------------
# Save RData
# --------------------------------------------------

save(
  net,
  MEs,
  moduleColors,
  moduleTraitCor,
  moduleTraitPvalue,
  traitData,
  meta,
  softPower,
  top_var_genes,
  file = file.path(
    outdir,
    "WGCNA_Result.RData"
  )
)

cat("WGCNA finished.\n")
cat("Samples:", nrow(datExpr), "\n")
cat("Genes:", ncol(datExpr), "\n")
cat("Modules:", length(unique(moduleColors)), "\n")
cat("Soft power:", softPower, "\n")