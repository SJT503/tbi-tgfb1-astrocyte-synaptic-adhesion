#!/usr/bin/env Rscript
# Step 5 confirmatory analysis -- INSPECT available objects
suppressMessages({library(Seurat)})
PD <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/processed_data"

inspect <- function(path){
  cat("\n========== ", path, " ==========\n")
  if(!file.exists(path)){ cat("  MISSING\n"); return(invisible()) }
  o <- readRDS(path)
  cat("  class:", class(o)[1], "\n")
  if(inherits(o,"Seurat")){
    cat("  cells:", ncol(o), " genes:", nrow(o), "\n")
    cat("  assays:", paste(Assays(o),collapse=","), " default:", DefaultAssay(o), "\n")
    md <- colnames(o@meta.data)
    cat("  meta cols:", paste(md, collapse=", "), "\n")
    for(cand in c("cell_type","celltype","CellType","seurat_annotations","annotation","major_celltype","cluster_annotation")){
      if(cand %in% md){ cat("  --", cand, "levels:\n"); print(table(o@meta.data[[cand]])) }
    }
    for(cand in c("dataset","orig.ident","GSE","batch","group","condition","sample","Group","Condition")){
      if(cand %in% md){ cat("  --", cand, ":\n"); print(table(o@meta.data[[cand]])) }
    }
  } else {
    cat("  str:\n"); str(o, max.level=2)
  }
}

inspect(file.path(PD,"integrated/mouse_integrated_annotated.rds"))
inspect(file.path(PD,"human_validation/gse209552_astro_validation.rds"))
inspect(file.path(PD,"spatial/visium_tbi_3dpi.rds"))
inspect(file.path(PD,"communication/cellchat_like_comparison.rds"))
cat("\nDONE inspect\n")
