# Script 101: 人源跨物种验证 — GSE209552 人类 TBI 星形谷氨酸处理面板
# ═══════════════════════════════════════════════════════════════════
# Script 100 在中度 CCI 小鼠（GSE226072）确认：星形丢失 GLT-1/GLAST/GS
# 本脚本：查人类 TBI（GSE209552）星形是否也丢失同一套基因 → 跨物种一致性
#   方法与 Script 100 一致：marker 定星形、深度降采样拉平、TBI vs CTL
#   人类基因符号大写：GLUL / SLC1A2(GLT-1) / SLC1A3(GLAST)

suppressPackageStartupMessages({ library(Seurat); library(dplyr) })

BASE   <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project"
INRDS  <- file.path(BASE, "wet_lab_data/results/human_validation_v2/processed.rds")
FIGOUT <- file.path(BASE, "wet_lab_data/results/figures/Fig_Homeostatic_Deepdive")
dir.create(FIGOUT, recursive=TRUE, showWarnings=FALSE)

PANEL <- list(
  "谷氨酸处理(↓预期)" = c("GLUL","SLC1A2","SLC1A3"),
  "反应性(↑预期)"     = c("GFAP","VIM","CD44","SERPINA3","OSMR","CP")
)
ALLG <- unlist(PANEL, use.names=FALSE)
MARKERS <- list(
  Neuron=c("SNAP25","SYT1","RBFOX3"), Astro=c("AQP4","GFAP","SLC1A2"),
  Oligo=c("PLP1","MBP","MOBP"), Micro=c("CSF1R","P2RY12","CX3CR1"),
  OPC=c("PDGFRA","CSPG4"), Endo=c("CLDN5","FLT1"))

cat("\n===== 加载人类 GSE209552 =====\n")
seu <- readRDS(INRDS)
seu$Group <- ifelse(grepl("HuBrainCTL", seu$sample_id), "CTL",
                    ifelse(grepl("MJ_TBI", seu$sample_id), "TBI", NA))
cat("分组:\n"); print(table(seu$Group))

# ── marker 定星形簇 ──
expr <- GetAssayData(seu, layer="data", assay="RNA")
cls <- sort(unique(seu$seurat_clusters))
scm <- sapply(names(MARKERS), function(ct) {
  g <- intersect(MARKERS[[ct]], rownames(expr))
  sapply(cls, function(cl) mean(as.matrix(expr[g, colnames(seu)[seu$seurat_clusters==cl], drop=FALSE])))
})
astro_cl <- cls[colnames(scm)[apply(scm,1,which.max)]=="Astro"]
cat(sprintf("星形簇: %s\n", paste(astro_cl, collapse=",")))
astro <- seu[, seu$seurat_clusters %in% astro_cl]
ae0 <- GetAssayData(astro, layer="data", assay="RNA")
pur <- mean(ae0["AQP4",]) / max(mean(ae0["SNAP25",]),1e-9)
cat(sprintf("人类星形 %d (CTL %d / TBI %d) | 纯度 AQP4/SNAP25=%.1f\n", ncol(astro),
            sum(astro$Group=="CTL"), sum(astro$Group=="TBI"), pur))

# ── 深度降采样拉平 ──
counts <- GetAssayData(astro, layer="counts", assay="RNA")
g <- astro$Group
med <- c(median(colSums(counts[,colnames(astro)[g=="TBI"],drop=FALSE])),
         median(colSums(counts[,colnames(astro)[g=="CTL"],drop=FALSE])))
target <- round(min(med)*0.9)
set.seed(42)
ds <- SampleUMI(counts, max.umi=target, upsample=FALSE, verbose=FALSE)
rownames(ds) <- rownames(counts); colnames(ds) <- colnames(counts)
dseu <- CreateSeuratObject(counts=ds)
dseu$Group <- g; dseu$sample_id <- astro$sample_id
dseu <- NormalizeData(dseu, verbose=FALSE)
de <- GetAssayData(dseu, layer="data", assay="RNA")

# ── 面板对比：per-cell + per-sample pseudobulk ──
cat("\n===== 人类星形谷氨酸/反应性面板（TBI vs CTL，深度匹配后）=====\n")
res <- lapply(ALLG, function(gene) {
  if (!gene %in% rownames(de)) return(data.frame(gene=gene, CTL=NA, TBI=NA, pct=NA, p_cell=NA, p_pb=NA))
  v_t <- as.numeric(de[gene, g=="TBI"]); v_c <- as.numeric(de[gene, g=="CTL"])
  p_cell <- tryCatch(wilcox.test(v_t, v_c)$p.value, error=function(e) NA)
  # per-sample pseudobulk
  pb <- data.frame(expr=as.numeric(de[gene,]), samp=dseu$sample_id, grp=g) %>%
    group_by(samp, grp) %>% summarise(m=mean(expr), .groups="drop")
  p_pb <- tryCatch(wilcox.test(pb$m[pb$grp=="TBI"], pb$m[pb$grp=="CTL"])$p.value,
                   error=function(e) NA)
  data.frame(gene=gene, CTL=round(mean(v_c),3), TBI=round(mean(v_t),3),
             pct=round((mean(v_t)-mean(v_c))/(mean(v_c)+1e-9)*100,1),
             p_cell=signif(p_cell,3), p_pb=round(p_pb,3))
})
tab <- bind_rows(res)
tab$category <- rep(names(PANEL), sapply(PANEL,length))
tab <- tab[,c("gene","category","CTL","TBI","pct","p_cell","p_pb")]
print(tab, row.names=FALSE)
write.csv(tab, file.path(FIGOUT,"human_glutamate_panel.csv"), row.names=FALSE)

cat("\n===== 跨物种判读 =====\n")
glu <- tab[tab$category=="谷氨酸处理(↓预期)",]
cat(sprintf("谷氨酸处理基因 %d 个：方向↓ %d/%d；per-sample 显著(p<0.05) %d/%d\n",
            nrow(glu), sum(glu$pct<0,na.rm=TRUE), nrow(glu),
            sum(glu$p_pb<0.05,na.rm=TRUE), nrow(glu)))
cat("\n小鼠中度 CCI（GSE226072）：GLUL/SLC1A2/SLC1A3 全部↓且显著\n")
cat("人类若同向↓ → GLT-1 丢失机制跨物种保守\n")
cat("\n========= 完成 =========\n")
