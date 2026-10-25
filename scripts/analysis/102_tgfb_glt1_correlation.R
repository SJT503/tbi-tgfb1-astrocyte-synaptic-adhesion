# Script 102: TGFβ 通路活性 vs GLT-1/GLAST/GS 相关性分析（计算地基）
# ═══════════════════════════════════════════════════════════════════
# 不是虚拟敲除（受体≠TF，工具不对）——是合法的相关性分析
# 问题：星形里 TGFβ 通路活性高的细胞，是不是 GLT-1/GLAST/GS 低？
#   负相关 → 支持"TGFβRI 信号压制谷氨酸处理"的方向（仍是相关，非因果）
# 数据：GSE226072（中度 CCI 小鼠）+ GSE209552（人类）星形
# 方法：Smad 靶基因签名评分 per 星形 → 与谷氨酸基因 Spearman 相关 + 高低分组对比

suppressPackageStartupMessages({ library(Seurat); library(dplyr) })

BASE <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project"
G226 <- file.path(BASE, "wet_lab_data/external_data/GSE226072")
HUM  <- file.path(BASE, "wet_lab_data/results/human_validation_v2/processed.rds")
FIGOUT <- file.path(BASE, "wet_lab_data/results/figures/Fig_Homeostatic_Deepdive")

# Smad/TGFβ 经典转录靶基因签名
SMAD_TARGETS_MM <- c("Serpine1","Smad7","Skil","Pmepa1","Junb","Id1","Id2","Id3",
                     "Tgfbi","Ltbp2","Cdkn1a","Smad6","Gadd45b","Bcl3")
GLU_MM <- c("Slc1a2","Slc1a3","Glul")
MARKERS_MM <- list(Neuron=c("Snap25","Syt1","Rbfox3"), Astro=c("Aqp4","Gfap","Slc1a2"),
  Oligo=c("Plp1","Mbp","Mobp"), Micro=c("Csf1r","P2ry12","Cx3cr1"),
  OPC=c("Pdgfra","Cspg4"), Endo=c("Cldn5","Flt1"))

corr_report <- function(astro, smad_genes, glu_genes, tgfbr1, label) {
  cat(sprintf("\n══════════ %s ══════════\n", label))
  expr <- GetAssayData(astro, layer="data", assay="RNA")
  sg <- intersect(smad_genes, rownames(expr))
  cat(sprintf("Smad 靶基因命中 %d/%d\n", length(sg), length(smad_genes)))
  astro <- AddModuleScore(astro, features=list(sg), name="TGFb_act", assay="RNA")
  tgfb_score <- astro$TGFb_act1
  br1 <- if (tgfbr1 %in% rownames(expr)) as.numeric(expr[tgfbr1,]) else rep(NA,ncol(astro))

  for (g in glu_genes) {
    if (!g %in% rownames(expr)) { cat(sprintf("  %s 不在\n",g)); next }
    gv <- as.numeric(expr[g,])
    # Smad 活性 vs 谷氨酸基因
    r1 <- suppressWarnings(cor.test(tgfb_score, gv, method="spearman"))
    # TGFBR1 表达 vs 谷氨酸基因
    r2 <- suppressWarnings(cor.test(br1, gv, method="spearman"))
    # TGFβ 高 vs 低（中位分）
    hi <- tgfb_score > median(tgfb_score)
    pct <- (mean(gv[hi]) - mean(gv[!hi])) / (mean(gv[!hi])+1e-9) * 100
    cat(sprintf("  %-7s | Smad活性 r=%+.3f(p=%.1e) | TGFBR1 r=%+.3f | TGFβ高组 vs 低组 %+.1f%%\n",
                g, r1$estimate, r1$p.value, r2$estimate, pct))
  }
}

# ═══════════ GSE226072 中度 CCI 小鼠 ═══════════
cat("\n########## GSE226072（中度 CCI 小鼠）##########\n")
S <- list(c("GSM7062619_D1_S1","Sham"),c("GSM7062620_D1_S2","CCI"),
          c("GSM7062623_D1_S5","Sham"),c("GSM7062624_D1_S6","CCI"),
          c("GSM7062627_D1_S9","Sham"),c("GSM7062628_D1_S10","CCI"),
          c("GSM7062631_D1_S13","Sham"),c("GSM7062632_D1_S14","CCI"))
objs <- lapply(S, function(s) {
  m <- ReadMtx(mtx=file.path(G226,paste0(s[1],"matrix.mtx.gz")),
               cells=file.path(G226,paste0(s[1],"barcodes.tsv.gz")),
               features=file.path(G226,paste0(s[1],"features.tsv.gz")), feature.column=2)
  o <- CreateSeuratObject(counts=m, min.cells=3, min.features=200); o$Group <- s[2]; o })
seu <- merge(objs[[1]], objs[-1]); seu <- JoinLayers(seu)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern="^mt-")
seu <- subset(seu, subset = nFeature_RNA>=200 & percent.mt<25)
seu <- NormalizeData(seu, verbose=FALSE); seu <- FindVariableFeatures(seu, nfeatures=2000, verbose=FALSE)
seu <- ScaleData(seu, verbose=FALSE); seu <- RunPCA(seu, npcs=30, verbose=FALSE)
seu <- FindNeighbors(seu, dims=1:30, verbose=FALSE); seu <- FindClusters(seu, resolution=0.5, verbose=FALSE)
ex <- GetAssayData(seu, layer="data", assay="RNA")
cls <- sort(unique(seu$seurat_clusters))
scm <- sapply(names(MARKERS_MM), function(ct) {
  g <- intersect(MARKERS_MM[[ct]], rownames(ex))
  sapply(cls, function(cl) mean(as.matrix(ex[g, colnames(seu)[seu$seurat_clusters==cl], drop=FALSE]))) })
astro_cl <- cls[colnames(scm)[apply(scm,1,which.max)]=="Astro"]
astro <- seu[, seu$seurat_clusters %in% astro_cl]
cat(sprintf("星形 %d\n", ncol(astro)))
corr_report(astro, SMAD_TARGETS_MM, GLU_MM, "Tgfbr1", "GSE226072 全部星形")
corr_report(astro[, astro$Group=="CCI"], SMAD_TARGETS_MM, GLU_MM, "Tgfbr1", "GSE226072 仅 CCI 星形")
rm(objs, seu, astro); gc(verbose=FALSE)

# ═══════════ GSE209552 人类 ═══════════
cat("\n########## GSE209552（人类 TBI）##########\n")
hum <- readRDS(HUM)
hum$Group <- ifelse(grepl("HuBrainCTL", hum$sample_id),"CTL",
             ifelse(grepl("MJ_TBI", hum$sample_id),"TBI",NA))
he <- GetAssayData(hum, layer="data", assay="RNA")
MARKERS_HS <- lapply(MARKERS_MM, toupper)
cls <- sort(unique(hum$seurat_clusters))
scm <- sapply(names(MARKERS_HS), function(ct) {
  g <- intersect(MARKERS_HS[[ct]], rownames(he))
  sapply(cls, function(cl) mean(as.matrix(he[g, colnames(hum)[hum$seurat_clusters==cl], drop=FALSE]))) })
astro_cl <- cls[colnames(scm)[apply(scm,1,which.max)]=="Astro"]
hastro <- hum[, hum$seurat_clusters %in% astro_cl]
cat(sprintf("人类星形 %d\n", ncol(hastro)))
corr_report(hastro, toupper(SMAD_TARGETS_MM), toupper(GLU_MM), "TGFBR1", "GSE209552 全部星形")
corr_report(hastro[, hastro$Group=="TBI"], toupper(SMAD_TARGETS_MM), toupper(GLU_MM), "TGFBR1", "GSE209552 仅 TBI 星形")

cat("\n判读：Smad活性 vs GLT-1/GLAST/GS 若负相关(r<0) + TGFβ高组比低组低 → 支持方向\n")
cat("\n========= 完成 =========\n")
