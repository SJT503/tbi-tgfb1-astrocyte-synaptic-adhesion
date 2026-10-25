#!/usr/bin/env Rscript
suppressMessages({library(Seurat)})
P <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project"
obj <- readRDS(file.path(P,"processed_data/integrated/mouse_integrated_annotated.rds"))
cat("class:", class(obj)[1], " dim:", paste(dim(obj),collapse=" x "), "\n")
cat("assays:", paste(Assays(obj),collapse=","), " default:", DefaultAssay(obj), "\n")
cat("--- meta columns ---\n"); print(colnames(obj@meta.data))
for (cc in colnames(obj@meta.data)) {
  v <- obj@meta.data[[cc]]
  if ((is.character(v)||is.factor(v)) && length(unique(v))<=30) { cat("[",cc,"]\n"); print(table(v)) }
}
cat("DONE\n")
