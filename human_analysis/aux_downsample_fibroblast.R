library(tidyverse)
library(Seurat)
library(DropletUtils)

here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
counts_dir <- file.path(here, "data", "fibroblast_counts")
res <- file.path(here, "results", "03_fibroblast_vs_SF", "downsampled")
dir.create(res, showWarnings = FALSE, recursive = TRUE)

m <- as(Matrix::readMM(file.path(counts_dir, "matrix.mtx")), "CsparseMatrix")
genes <- read_lines(file.path(counts_dir, "genes.tsv"))
meta.data <- read.delim(file.path(counts_dir, "cells.tsv"), row.names = 1)
dimnames(m) <- list(genes, rownames(meta.data))

roi <- meta.data$roi_name
med <- tapply(Matrix::colSums(m), roi, median)
prop <- min(med) / med
m.downsample <- do.call(cbind, lapply(names(prop), \(s) {
  set.seed(0)
  downsampleMatrix(m[, roi == s], prop = prop[[s]], bycol = FALSE)
}))[, rownames(meta.data)]
print(cbind(median = med, prop = prop, median_after = tapply(Matrix::colSums(m.downsample), roi, median)))

obj <- CreateSeuratObject(m.downsample, meta.data = meta.data[, c("roi_name", "group")])
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 100)

pairs <- list(c("Xray_Late", "SF1"), c("Xray_Late", "SF2"), c("Xray_Late", "SF"),
              c("Cion_Late", "SF1"), c("Cion_Late", "SF2"), c("Cion_Late", "SF"),
              c("Xray_Late", "Cion_Late"),
              c("Xray_Under", "SF1"), c("Xray_Under", "SF2"), c("Xray_Under", "SF"),
              c("Xray_Under", "Xray_Late"), c("Xray_Under", "Cion_Late"))

deg <- function(a, b) {
  key <- function(g) if (g == "SF") "group" else "roi_name"
  cells.a <- colnames(obj)[obj[[key(a)]][, 1] == a]
  cells.b <- colnames(obj)[obj[[key(b)]][, 1] == b]
  r <- FindMarkers(obj, ident.1 = cells.a, ident.2 = cells.b, test.use = "wilcox",
                   logfc.threshold = 0, min.pct = 0, return.thresh = 1)
  r$p_val_adj <- p.adjust(r$p_val, method = "BH")
  r$sig <- r$p_val_adj < 0.05 & abs(r$avg_log2FC) > 0.5 & pmax(r$pct.1, r$pct.2) >= 0.1
  r
}

results <- lapply(pairs, \(p) deg(p[1], p[2]))
names(results) <- sapply(pairs, \(p) paste(p[1], "vs", p[2]))
for (n in names(results)) write.csv(results[[n]], file.path(res, paste0("deg_", gsub(" vs ", "_vs_", n), ".csv")))

summary <- tibble(comparison = names(results),
                  up_in_A = sapply(results, \(r) sum(r$sig & r$avg_log2FC > 0)),
                  up_in_B = sapply(results, \(r) sum(r$sig & r$avg_log2FC < 0)))
write_csv(summary, file.path(res, "summary.csv"))
print(summary, n = Inf)

both <- function(x, direction) {
  s <- function(r) rownames(r)[r$sig & sign(r$avg_log2FC) == direction]
  sort(intersect(s(results[[paste(x, "vs SF1")]]), s(results[[paste(x, "vs SF2")]])))
}
for (x in c("Xray_Under", "Xray_Late", "Cion_Late")) {
  cat(x, "> SF1 and SF2:", both(x, 1), "\n")
  cat(x, "< SF1 and SF2:", both(x, -1), "\n")
}

fib.genes <- intersect(c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "FN1", "POSTN", "ELN", "FBLN5",
                         "MFAP5", "SERPINE1", "SNAI2", "CCN2", "CTHRC1", "LOX", "TNC", "THBS1", "TGFB1",
                         "TIMP1", "SPARC"), rownames(obj))
lfc <- sapply(results, \(r) paste0(round(r[fib.genes, "avg_log2FC"], 2), ifelse(r[fib.genes, "sig"], " *", "")))
rownames(lfc) <- fib.genes
write.csv(lfc, file.path(res, "fibrosis_genes_log2fc.csv"))
print(lfc, quote = FALSE)
