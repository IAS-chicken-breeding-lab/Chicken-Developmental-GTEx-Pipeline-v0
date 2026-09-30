```r
#!/usr/bin/env Rscript

library(edgeR)
library(data.table)

tis <- commandArgs(trailingOnly = TRUE)[1]

counts <- read.table(paste0(tis, ".counts"),
                     header = TRUE, row.names = 1)

tpm <- read.table(paste0(tis, ".tpm"),
                  header = TRUE, row.names = 1)

tpm <- tpm[rownames(counts), colnames(counts)]

dge <- calcNormFactors(DGEList(counts = counts), method = "TMM")
tmm <- cpm(dge, normalized.lib.sizes = TRUE)

keep <- rowSums(counts >= 6) >= 0.2 * ncol(counts) &
        rowSums(tpm >= 0.1) >= 0.2 * ncol(tpm)

tmm <- tmm[keep, ]

tmm_inv <- t(apply(tmm, 1, function(x)
  qnorm((rank(x) - 0.5) / length(x))
))

fwrite(
  data.frame(gene_id = rownames(tmm), tmm),
  paste0(tis, "_TMM_pass.txt"),
  sep = "\t"
)

bed <- fread("geno_eqtl.bed", header = FALSE)
setnames(bed, c("#Chr", "start", "end", "gene_id"))

bed <- bed[
  `#Chr` %in% as.character(1:39) &
  gene_id %in% rownames(tmm_inv)
]

bed_out <- cbind(
  bed,
  tmm_inv[bed$gene_id, , drop = FALSE]
)

fwrite(bed_out, paste0(tis, ".bed"), sep = "\t")

system(paste0("bgzip ", tis, ".bed"))
system(paste0("tabix -p bed ", tis, ".bed.gz"))
```
