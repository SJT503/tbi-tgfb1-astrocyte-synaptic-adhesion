#!/usr/bin/env Rscript
# 实查 severity_integrated_allcells.rds 是否含 CCI 皮层 7d 数据 (作 SFig 对照)
suppressMessages({library(Seurat)})
RDS <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/processed_data/v2/severity_integrated_allcells.rds"
cat("probe:", RDS, "\n")
seu <- readRDS(RDS)
cat("class:", paste(class(seu), collapse=","), "\n")
cat("cells:", ncol(seu), "  features:", nrow(seu), "\n")
mt <- seu@meta.data
cat("meta cols:", paste(colnames(mt), collapse=","), "\n\n")
for (cl in c("dataset","severity","condition","time_point","timepoint","sample","celltype","celltype_auto","stim","group")) {
  if (cl %in% colnames(mt)) {
    v <- mt[[cl]]; u <- unique(v)
    if (length(u) <= 50) { cat(sprintf("[%s] (%d):\n", cl, length(u))); print(table(v)); cat("\n") }
  }
}
cat("\n=== celltype-related columns ===\n")
ct_cols <- grep("celltype|cell_type|annotation", colnames(mt), value=TRUE, ignore.case=TRUE)
cat("found:", paste(ct_cols, collapse=","), "\n")
if (length(ct_cols)==0) {
  cat("NO celltype column — need to derive from clusters or use markers\n")
  cat("\nseurat_clusters distribution:\n")
  if ("seurat_clusters" %in% colnames(mt)) print(table(mt$seurat_clusters))
}
cat("\n=== dataset x time_point cross-tab ===\n")
tcol <- if ("time_point" %in% colnames(mt)) "time_point" else if ("timepoint" %in% colnames(mt)) "timepoint" else NULL
if (!is.null(tcol)) print(table(mt$dataset, mt[[tcol]]))
