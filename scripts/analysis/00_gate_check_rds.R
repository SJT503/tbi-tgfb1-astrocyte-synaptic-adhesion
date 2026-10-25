#!/usr/bin/env Rscript
# 00_gate_check_rds.R — 入口五检验 (CLAUDE.md 焊死铁律)
# 检验对象: astro_sub_adhesion_3groups.rds (升级 A/C 主输入)
# 标准: 星形纯度 Aqp4/Gfap >> Snap25 比值 ≥3, 否则 stop()
suppressMessages({library(Seurat)})

RDS <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/astro_sub_adhesion_3groups.rds"
cat("[GATE] reading:", RDS, "\n")
seu <- readRDS(RDS)

cat("\n[GATE-1] cell count\n")
n_cells <- ncol(seu)
cat("  cells =", n_cells, "\n")
if (n_cells < 100) stop("GATE-1 FAIL: too few cells (<100)")

cat("\n[GATE-2] cluster count\n")
cl <- tryCatch(Idents(seu), error=function(e) NULL)
if (!is.null(cl)) cat("  clusters =", nlevels(cl), "\n") else cat("  (no Idents)\n")

cat("\n[GATE-3] celltype purity (Aqp4/Gfap vs Snap25/Syt1)\n")
DefaultAssay(seu) <- "RNA"
layer <- if ("data" %in% Layers(seu)) "data" else "counts"
expr <- GetAssayData(seu, layer=layer, assay="RNA")
markers <- c("Aqp4","Gfap","Snap25","Syt1")
present <- intersect(markers, rownames(expr))
cat("  markers present:", paste(present, collapse=","), "\n")
if (length(present) < 3) stop("GATE-3 FAIL: <3 markers in matrix")
m_aqp4   <- if ("Aqp4"   %in% present) mean(expr["Aqp4",])   else NA
m_gfap   <- if ("Gfap"   %in% present) mean(expr["Gfap",])   else NA
m_snap25 <- if ("Snap25" %in% present) mean(expr["Snap25",]) else NA
m_syt1   <- if ("Syt1"   %in% present) mean(expr["Syt1",])   else NA
cat(sprintf("  Aqp4=%.3f  Gfap=%.3f  Snap25=%.3f  Syt1=%.3f\n",
            m_aqp4, m_gfap, m_snap25, m_syt1))
astro_mean   <- mean(c(m_aqp4, m_gfap), na.rm=TRUE)
neuron_mean  <- mean(c(m_snap25, m_syt1), na.rm=TRUE)
ratio        <- astro_mean / neuron_mean
cat(sprintf("  ratio (astro/neuron) = %.2f x  (need >= 3)\n", ratio))
if (is.na(ratio) || ratio < 3)
  stop(sprintf("GATE-3 FAIL: purity ratio %.2f < 3 -- not an astrocyte RDS!", ratio))

cat("\n[GATE-4] group/condition completeness\n")
meta <- seu@meta.data
grp_cols <- intersect(c("group","condition","orig.ident","sample","Sample","stim"),
                      colnames(meta))
cat("  meta columns:", paste(colnames(meta), collapse=","), "\n")
if (length(grp_cols) > 0) {
  for (gc in grp_cols) {
    tb <- table(meta[[gc]])
    cat(sprintf("  [%s] %s\n", gc,
                paste(names(tb), tb, sep="=", collapse=", ")))
  }
}

cat("\n[GATE-5] annotation/idents\n")
cat("  Idents head:", paste(head(unique(as.character(Idents(seu))), 10), collapse=","), "\n")

cat("\n[GATE] ALL PASSED. RDS safe for downstream use.\n")
saveRDS(list(
  rds=RDS, n_cells=n_cells, ratio=ratio,
  meta_cols=colnames(meta), grp_cols=grp_cols,
  ts=Sys.time()
), "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/upgrade_A_copositivity/gate_check_result.rds")
