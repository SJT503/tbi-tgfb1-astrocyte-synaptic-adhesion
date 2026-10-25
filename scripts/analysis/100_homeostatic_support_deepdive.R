# Script 100: "星形稳态支持丢失"机制深挖 — 中度 CCI 是否成立 + KO 能否恢复
# ═══════════════════════════════════════════════════════════════════
# Script 99 发现 robust 跨数据集星形 TBI 签名 = 反应性↑ + 稳态支持↓
#   但 robust 签名主要来自轻度 Drop-seq 数据集
# 本脚本：在两个关键数据集里直接查这套基因面板的方向
#   ① GSE226072（中度 CCI snRNA）— Sham vs CCI：中度损伤里成不成立？
#   ② 自测深度匹配数据 — Sham vs CCI（损伤）+ CCI vs KO（救治）：KO 能否恢复？
# 不卡 DEG 阈值，直接看面板基因的方向 + 显著性
# ⚠️ 自测 KO 对比仍受批次混杂影响（深度已匹配但批次没修）→ 提示性，非定论

suppressPackageStartupMessages({ library(Seurat); library(dplyr) })

BASE   <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project"
G226   <- file.path(BASE, "wet_lab_data/external_data/GSE226072")
DM     <- file.path(BASE, "wet_lab_data/results/seurat/depth_matched/astrocyte_depthmatched.rds")
FIGOUT <- file.path(BASE, "wet_lab_data/results/figures/Fig_Homeostatic_Deepdive")
dir.create(FIGOUT, recursive=TRUE, showWarnings=FALSE)

# ── 基因面板（Script 99 robust 信号 + 经典星形稳态/反应性）─────────────────────
PANEL <- list(
  "反应性(↑预期)"   = c("Gfap","Vim","Cd44","Lcn2","Serpina3n","Osmr","Cp","Hbegf"),
  "谷氨酸处理(↓预期)" = c("Glul","Slc1a2","Slc1a3"),
  "缝隙连接(↓预期)" = c("Gja1","Gjb6"),
  "代谢水通道(↓预期)" = c("Aldoc","Aqp4","Mlc1","Ldhb"),
  "其他稳态(↓预期)" = c("Cldn10","Sparcl1","Fam107a","Tspan7")
)
ALLG <- unlist(PANEL, use.names=FALSE)
MARKERS <- list(
  Neuron=c("Snap25","Syt1","Rbfox3"), Astro=c("Aqp4","Gfap","Slc1a2"),
  Oligo=c("Plp1","Mbp","Mobp"), Micro=c("Csf1r","P2ry12","Cx3cr1"),
  OPC=c("Pdgfra","Cspg4"), Endo=c("Cldn5","Flt1"))

panel_compare <- function(expr, g1cells, g2cells, g1name, g2name) {
  # g1 = 损伤/KO, g2 = 对照/CCI；返回每个面板基因的方向 + p
  res <- lapply(ALLG, function(g) {
    if (!g %in% rownames(expr)) return(data.frame(gene=g, m1=NA, m2=NA, pct=NA, p=NA))
    v1 <- as.numeric(expr[g, g1cells]); v2 <- as.numeric(expr[g, g2cells])
    p <- tryCatch(wilcox.test(v1, v2)$p.value, error=function(e) NA)
    data.frame(gene=g, m1=mean(v1), m2=mean(v2),
               pct=round((mean(v1)-mean(v2))/(mean(v2)+1e-9)*100,1), p=p)
  })
  d <- bind_rows(res)
  names(d) <- c("gene", g1name, g2name, paste0(g1name,"_vs_",g2name,"_pct"),
                paste0(g1name,"_vs_",g2name,"_p"))
  d
}

# ═══════════ ① GSE226072（中度 CCI snRNA）═══════════
cat("\n========== GSE226072 中度 CCI — 星形稳态面板 ==========\n")
S <- list(c("GSM7062619_D1_S1","Sham"),c("GSM7062620_D1_S2","CCI"),
          c("GSM7062623_D1_S5","Sham"),c("GSM7062624_D1_S6","CCI"),
          c("GSM7062627_D1_S9","Sham"),c("GSM7062628_D1_S10","CCI"),
          c("GSM7062631_D1_S13","Sham"),c("GSM7062632_D1_S14","CCI"))
objs <- lapply(S, function(s) {
  m <- ReadMtx(mtx=file.path(G226,paste0(s[1],"matrix.mtx.gz")),
               cells=file.path(G226,paste0(s[1],"barcodes.tsv.gz")),
               features=file.path(G226,paste0(s[1],"features.tsv.gz")), feature.column=2)
  o <- CreateSeuratObject(counts=m, min.cells=3, min.features=200)
  o$Group <- s[2]; o$samp <- s[1]; o
})
seu <- merge(objs[[1]], objs[-1]); seu <- JoinLayers(seu)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern="^mt-")
seu <- subset(seu, subset = nFeature_RNA>=200 & percent.mt<25)
seu <- NormalizeData(seu, verbose=FALSE)
seu <- FindVariableFeatures(seu, nfeatures=2000, verbose=FALSE)
seu <- ScaleData(seu, verbose=FALSE); seu <- RunPCA(seu, npcs=30, verbose=FALSE)
seu <- FindNeighbors(seu, dims=1:30, verbose=FALSE)
seu <- FindClusters(seu, resolution=0.5, verbose=FALSE)
expr <- GetAssayData(seu, layer="data", assay="RNA")
cls <- sort(unique(seu$seurat_clusters))
scm <- sapply(names(MARKERS), function(ct) {
  g <- intersect(MARKERS[[ct]], rownames(expr))
  sapply(cls, function(cl) mean(as.matrix(expr[g, colnames(seu)[seu$seurat_clusters==cl], drop=FALSE])))
})
astro_cl <- cls[colnames(scm)[apply(scm,1,which.max)]=="Astro"]
astro <- seu[, seu$seurat_clusters %in% astro_cl]
ae <- GetAssayData(astro, layer="data", assay="RNA")
cat(sprintf("GSE226072 星形 %d (Sham %d / CCI %d)\n", ncol(astro),
            sum(astro$Group=="Sham"), sum(astro$Group=="CCI")))
g226 <- panel_compare(ae, colnames(astro)[astro$Group=="CCI"],
                      colnames(astro)[astro$Group=="Sham"], "G226_CCI", "G226_Sham")
rm(objs, seu, astro); gc(verbose=FALSE)

# ═══════════ ② 自测深度匹配星形 ═══════════
cat("\n========== 自测深度匹配星形 — 损伤 + 救治 ==========\n")
sf <- readRDS(DM)
se <- GetAssayData(sf, layer="data", assay="RNA")
g <- sf$stim
cat(sprintf("自测星形 %d (Sham %d / CCI_Scr %d / KO %d)\n", ncol(sf),
            sum(g=="Sham"), sum(g=="CCI_Scr"), sum(g=="CCI_TGFbetaR1_KO")))
self_inj <- panel_compare(se, colnames(sf)[g=="CCI_Scr"], colnames(sf)[g=="Sham"],
                          "self_CCI", "self_Sham")
self_ko  <- panel_compare(se, colnames(sf)[g=="CCI_TGFbetaR1_KO"], colnames(sf)[g=="CCI_Scr"],
                          "self_KO", "self_CCI2")

# ═══════════ 汇总 ═══════════
cat("\n========== 汇总：星形稳态支持面板 ==========\n")
panel_cat <- data.frame(gene=ALLG,
                        category=rep(names(PANEL), sapply(PANEL,length)))
merged <- panel_cat %>%
  left_join(g226[,c("gene","G226_CCI_vs_G226_Sham_pct","G226_CCI_vs_G226_Sham_p")], by="gene") %>%
  left_join(self_inj[,c("gene","self_CCI_vs_self_Sham_pct","self_CCI_vs_self_Sham_p")], by="gene") %>%
  left_join(self_ko[,c("gene","self_KO_vs_self_CCI2_pct","self_KO_vs_self_CCI2_p")], by="gene")
names(merged) <- c("gene","category","G226_中度CCI_%","G226_p",
                   "自测_损伤_%","自测损伤_p","自测_KO救治_%","自测KO_p")
merged <- merged %>% mutate(across(where(is.numeric), ~round(.,3)))
print(merged, row.names=FALSE)
write.csv(merged, file.path(FIGOUT,"homeostatic_panel_3comparisons.csv"), row.names=FALSE)

# ═══════════ 一致性判读 ═══════════
cat("\n========== 判读 ==========\n")
chk <- function(sub, pct_col, expect) {
  v <- merged[[pct_col]][merged$category==sub]
  v <- v[!is.na(v)]
  if (expect=="up") sum(v>0) else sum(v<0)
}
for (cat_name in names(PANEL)) {
  exp_dir <- if (grepl("↑", cat_name)) "up" else "down"
  n_tot <- sum(merged$category==cat_name & !is.na(merged$`G226_中度CCI_%`))
  cat(sprintf("\n[%s] 共 %d 基因，预期方向 %s\n", cat_name, n_tot, exp_dir))
  cat(sprintf("  GSE226072 中度CCI 符合预期: %d/%d\n", chk(cat_name,"G226_中度CCI_%",exp_dir), n_tot))
  cat(sprintf("  自测 损伤      符合预期: %d/%d\n", chk(cat_name,"自测_损伤_%",exp_dir), n_tot))
  # KO 救治：应与损伤反向
  ko_dir <- if (exp_dir=="up") "down" else "up"
  cat(sprintf("  自测 KO 救治(应反向 %s): %d/%d\n", ko_dir, chk(cat_name,"自测_KO救治_%",ko_dir), n_tot))
}
cat("\n========= 完成 =========\n")
