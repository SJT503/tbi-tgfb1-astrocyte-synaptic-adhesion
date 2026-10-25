#!/usr/bin/env Rscript
# C6 probe — 看 visium_tbi_3dpi.rds 结构, 判断空间通讯加强方案是否可行
suppressMessages({library(Seurat)})
RDS <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/processed_data/spatial/visium_tbi_3dpi.rds"
x <- readRDS(RDS)
cat("class:", paste(class(x), collapse=","), "\n")
if (inherits(x, "Seurat")) {
  cat("spots:", ncol(x), " features:", nrow(x), "\n")
  cat("assays:", paste(Assays(x), collapse=","), "\n")
  cat("default:", DefaultAssay(x), "\n")
  cat("images:", paste(Images(x), collapse=","), "\n")
  mt <- x@meta.data
  cat("meta cols:", paste(colnames(mt), collapse=","), "\n")
  for (cl in colnames(mt)) {
    v <- mt[[cl]]
    if (is.character(v) || is.factor(v)) {
      u <- unique(v); if (length(u) <= 20) { cat(sprintf("[%s] (%d):\n", cl, length(u))); print(table(v)) }
    }
  }
  cat("\nIdents head:", paste(head(unique(as.character(Idents(x))), 15), collapse=","), "\n")
}
