#!/usr/bin/env Rscript
suppressMessages({library(Seurat)})
RDS <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/astro_sub_adhesion_3groups.rds"
x <- readRDS(RDS)
cat("class:", paste(class(x), collapse=","), "\n")
cat("length:", length(x), "\n")
if (is.list(x)) {
  cat("names:", paste(names(x), collapse=" | "), "\n")
  for (nm in names(x)) {
    el <- x[[nm]]
    cat(sprintf("  %s -> %s (len %d)\n", nm,
                paste(class(el), collapse=","), length(el)))
    if (inherits(el, "Seurat")) {
      cat(sprintf("    cells=%d, features=%d, assays=%s, default=%s\n",
                  ncol(el), nrow(el),
                  paste(Assays(el), collapse=","),
                  DefaultAssay(el)))
      mt <- el@meta.data
      cat("    meta cols:", paste(colnames(mt), collapse=","), "\n")
    }
  }
} else if (inherits(x, "Seurat")) {
  cat("Seurat: cells=", ncol(x), " features=", nrow(x), "\n", sep="")
  cat("meta cols:", paste(colnames(x@meta.data), collapse=","), "\n")
}
