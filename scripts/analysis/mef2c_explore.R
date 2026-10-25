#!/usr/bin/env Rscript
# MEF2C 验证 - 第0步: 探索公共数据对象结构
suppressMessages({library(Seurat)})
P <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project"

inspect <- function(path, sp_markers) {
  cat("\n========== ", basename(path), " ==========\n")
  if (!file.exists(path)) { cat("  FILE NOT FOUND\n"); return(invisible()) }
  obj <- readRDS(path)
  cat("class:", class(obj)[1], "\n")
  if (inherits(obj, "Seurat")) {
    cat("dim (genes x cells):", paste(dim(obj), collapse=" x "), "\n")
    cat("assays:", paste(Assays(obj), collapse=", "), "\n")
    da <- DefaultAssay(obj)
    cat("DefaultAssay:", da, "\n")
    cat("layers:", paste(tryCatch(Layers(obj), error=function(e) "NA"), collapse=", "), "\n")
    cat("--- metadata columns ---\n")
    print(colnames(obj@meta.data))
    for (cc in colnames(obj@meta.data)) {
      v <- obj@meta.data[[cc]]
      if ((is.character(v) || is.factor(v)) && length(unique(v)) <= 20) {
        cat("  [", cc, "]: "); print(table(v))
      }
    }
    # purity
    expr <- tryCatch(GetAssayData(obj, assay="RNA", layer="data"),
                     error=function(e) GetAssayData(obj, assay=da, layer="data"))
    cat("--- purity / gene presence ---\n")
    for (g in sp_markers) {
      hit <- rownames(expr)[toupper(rownames(expr)) == toupper(g)]
      if (length(hit)>0) cat("  ", g, "(", hit[1], ") mean(data):", round(mean(expr[hit[1],]),4), "\n")
      else cat("  ", g, ": NOT in matrix\n")
    }
    # MEF2C + adhesion genes presence
    cat("--- target genes presence ---\n")
    for (g in c("Mef2c","MEF2C","Nrxn1","NRXN1","Nrxn3","Cadm1","CADM1","Cntn1","CNTN1","Cdh2","CDH2","Tgfbr1","TGFBR1")) {
      hit <- rownames(expr)[rownames(expr) == g]
      if (length(hit)>0) cat("  ", g, ": present, mean(data)=", round(mean(expr[g,]),4), "\n")
    }
  } else {
    cat("str:\n"); print(utils::head(obj))
  }
  invisible()
}

inspect(file.path(P,"processed_data/human_validation/gse226207_astro_validation.rds"),
        c("Aqp4","Gfap","Snap25"))
inspect(file.path(P,"processed_data/human_validation/gse209552_astro_validation.rds"),
        c("AQP4","GFAP","SNAP25"))

cat("\n========== regulons_pruned.rds ==========\n")
rp <- file.path(P,"scenic_output/regulons_pruned.rds")
if (file.exists(rp)) {
  reg <- readRDS(rp)
  cat("class:", class(reg)[1], " length/dim:", paste(dim(as.matrix(reg)), collapse="x"), "\n")
  nm <- tryCatch(names(reg), error=function(e) NULL)
  if (!is.null(nm)) { cat("n regulons:", length(nm), "\n")
    cat("Mef2c regulon present:", any(grepl("Mef2c", nm, ignore.case=TRUE)), "\n")
    cat("first 30:\n"); print(utils::head(nm,30)) }
}
cat("\n========== aucell_gse226207.rds ==========\n")
ag <- file.path(P,"scenic_output/aucell_gse226207.rds")
if (file.exists(ag)) {
  au <- readRDS(ag)
  cat("class:", class(au)[1], "\n")
  m <- tryCatch(as.matrix(au), error=function(e) NULL)
  if (!is.null(m)) { cat("dim:", paste(dim(m), collapse=" x "), "\n")
    rn <- rownames(m); cat("Mef2c regulon row:", any(grepl("Mef2c", rn, ignore.case=TRUE)), "\n")
    cat("first 25 regulon names:\n"); print(utils::head(rn,25)) }
}
cat("\nDONE explore\n")
