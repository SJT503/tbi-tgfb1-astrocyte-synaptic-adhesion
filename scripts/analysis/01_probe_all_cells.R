#!/usr/bin/env Rscript
# probe wet_lab_data/results/seurat/all_cells_annotated.rds
# 目的: 确认是自测三组 Sham/CCI_Scr/CCI_KO 全细胞 Seurat, 含星形+神经元
suppressMessages({library(Seurat)})
RDS <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/wet_lab_data/results/seurat/all_cells_annotated.rds"
cat("probe:", RDS, "\n")
seu <- readRDS(RDS)
cat("class:", paste(class(seu), collapse=","), "\n")
if (!inherits(seu, "Seurat")) { cat("NOT Seurat\n"); quit() }
cat("cells:", ncol(seu), " features:", nrow(seu), "\n")
cat("assays:", paste(Assays(seu), collapse=","), " default:", DefaultAssay(seu), "\n")
mt <- seu@meta.data
cat("\nmeta cols:", paste(colnames(mt), collapse=","), "\n")
for (cl in colnames(mt)) {
  v <- mt[[cl]]
  if (is.character(v) || is.factor(v)) {
    u <- unique(v); nu <- length(u)
    if (nu <= 30) {
      tb <- table(v)
      cat(sprintf("  [%s] (%d levels): %s\n", cl, nu,
                  paste(names(tb), tb, sep="=", collapse=", ")))
    } else {
      cat(sprintf("  [%s] (%d levels, sample): %s ...\n", cl, nu,
                  paste(head(u, 5), collapse=", ")))
    }
  }
}
cat("\nIdents head:", paste(head(unique(as.character(Idents(seu))), 20), collapse=","), "\n")
