#!/usr/bin/env Rscript
# Probe mouse_integrated_annotated.rds — meta 结构, 时间点分布, neuron 子集
suppressMessages({library(Seurat)})
RDS <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/processed_data/integrated/mouse_integrated_annotated.rds"
cat("probe:", RDS, "\n")
seu <- readRDS(RDS)
cat("class:", paste(class(seu), collapse=","), "\n")
cat("cells:", ncol(seu), "  features:", nrow(seu), "\n")
cat("assays:", paste(Assays(seu), collapse=","), "  default:", DefaultAssay(seu), "\n\n")

mt <- seu@meta.data
cat("meta cols:", paste(colnames(mt), collapse=", "), "\n\n")

for (cl in colnames(mt)) {
  v <- mt[[cl]]
  if (is.character(v) || is.factor(v)) {
    u <- unique(v); nu <- length(u)
    if (nu <= 20) {
      tb <- table(v)
      cat(sprintf("[%s] (%d levels):\n", cl, nu))
      print(tb); cat("\n")
    }
  }
}

# cross-tab: dataset x condition x time (if any)
cands <- intersect(c("dataset","orig.ident","condition","time","timepoint","stim","group","celltype","cell_type"),
                   colnames(mt))
cat("===== cross-tabs =====\n")
if ("celltype" %in% colnames(mt) || "cell_type" %in% colnames(mt)) {
  ctcol <- if ("celltype" %in% colnames(mt)) "celltype" else "cell_type"
  cat("celltype distribution:\n"); print(table(mt[[ctcol]])); cat("\n")
}
for (a in cands) for (b in cands) if (a < b) {
  cat(sprintf("--- %s x %s ---\n", a, b))
  print(table(mt[[a]], mt[[b]])); cat("\n")
}
