#!/usr/bin/env Rscript
# 人源 GSE209552 神经元 anoikis signature 分析
# 16 个样本 (11 TBI + 5 CTL) 缓存 RDS, 每个含全细胞类型
suppressMessages({library(Seurat); library(Matrix)})
set.seed(1)

CACHE <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/processed_data/human_validation/gse209552_cache"
OUT   <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/anoikis_out"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

rds_files <- list.files(CACHE, pattern="\\.rds$", full.names=TRUE)
cat("缓存样本:", length(rds_files), "\n")

# 分组: MJ_TBI_ 开头 = TBI; TBI_HuBrainCTL_ 开头 = CTL
get_group <- function(f){
  b <- basename(f)
  if (grepl("HuBrainCTL", b)) "CTL" else "TBI"
}
get_sample <- function(f) sub("\\.rds$","",basename(f))

# 人源标志物
M <- list(
  neuron  = c("SNAP25","SYT1","RBFOX3","MEG3","SYN1","STMN2"),
  astro   = c("AQP4","GFAP","SLC1A3","SLC1A2","ALDH1L1","GJA1"),
  micro   = c("CX3CR1","HEXB","C1QA","CSF1R","P2RY12"),
  oligo   = c("PLP1","MBP","MOG","MOBP"),
  opc     = c("PDGFRA","CSPG4"),
  vasc    = c("CLDN5","FLT1","PDGFRB"))

# anoikis (人源大写 symbol)
survival_genes  <- c("PTK2","AKT1","AKT2","AKT3","BCL2","MCL1","BCL2L1","BIRC5","MAPK1")
death_genes     <- c("BCL2L11","BMF","BAX","BAD","CASP9","CASP3","CYCS","APAF1")

objs <- list()
for (f in rds_files){
  o <- readRDS(f)
  if (!inherits(o, "Seurat")) { cat("  skip non-Seurat:", basename(f), "\n"); next }
  o$group  <- get_group(f)
  o$sample <- get_sample(f)
  # 统一 RNA assay
  if (!"RNA" %in% Assays(o)) next
  DefaultAssay(o) <- "RNA"
  objs[[get_sample(f)]] <- o
}
cat("有效 Seurat 对象:", length(objs), "\n")
obj <- merge(objs[[1]], objs[-1])
obj <- JoinLayers(obj, assay="RNA")
cat("合并细胞:", ncol(obj), "  分组:\n")
print(table(obj$group))

# 标准化(如未做)
if (max(GetAssayData(obj, layer="data")[,1:50]) == 0)
  obj <- NormalizeData(obj, verbose=FALSE)

# 细胞分类
sc <- sapply(names(M), function(k){
  g <- intersect(M[[k]], rownames(obj))
  if (length(g)==0) return(rep(0, ncol(obj)))
  Matrix::colMeans(GetAssayData(obj, layer="data")[g,,drop=FALSE])
})
obj$celltype <- colnames(sc)[max.col(sc, ties.method="first")]
cat("\n细胞类型分布:\n"); print(table(obj$celltype))

# 神经元子集 + 纯度门禁
neu <- subset(obj, subset = celltype=="neuron")
expr <- GetAssayData(neu, layer="data")
pr <- mean(expr["SNAP25",])/max(mean(expr["AQP4",]),1e-6)
cat(sprintf("\n神经元子集: n=%d  纯度 SNAP25/AQP4=%.1f\n", ncol(neu), pr))
if (pr < 3) stop("PURITY FAIL: human neuron SNAP25 not >> AQP4")

# 计算 anoikis signature score
neu <- AddModuleScore(neu, features=list(survival_genes), name="survival_score", search=TRUE)
neu <- AddModuleScore(neu, features=list(death_genes),    name="death_score",    search=TRUE)
neu$survival <- neu$survival_score1
neu$death    <- neu$death_score1
neu$anoikis  <- neu$death - neu$survival

# 每细胞 + 每样本汇总
cell_df <- data.frame(sample=neu$sample, group=neu$group,
                      survival=neu$survival, death=neu$death, anoikis=neu$anoikis)
write.csv(cell_df, file.path(OUT,"anoikis_human_per_cell.csv"), row.names=FALSE)

agg <- aggregate(cbind(survival,death,anoikis) ~ sample+group, cell_df, mean)
agg$group <- factor(agg$group, levels=c("CTL","TBI"))
write.csv(agg, file.path(OUT,"anoikis_human_per_sample.csv"), row.names=FALSE)

cat("\n=== 每样本汇总 ===\n"); print(agg[order(agg$group),], row.names=FALSE)

cat("\n=== 各组均值±SE ===\n")
for (s in c("survival","death","anoikis")){
  m  <- tapply(agg[[s]], agg$group, mean)
  se <- tapply(agg[[s]], agg$group, function(x) sd(x)/sqrt(length(x)))
  cat(s,":\n"); for (g in names(m)) cat(sprintf("  %-4s  %+.3f ± %.3f (n=%d)\n",g,m[g],se[g],sum(agg$group==g)))
  v1 <- agg[[s]][agg$group=="CTL"]; v2 <- agg[[s]][agg$group=="TBI"]
  p <- wilcox.test(v1, v2)$p.value
  cat(sprintf("    CTL vs TBI: Mann-Whitney p=%.4f\n", p))
}

# 关键基因表达 (用于 heatmap)
key_genes <- c(death_genes, survival_genes)
key_genes <- intersect(key_genes, rownames(neu))
gene_df <- data.frame(t(as.matrix(GetAssayData(neu, layer='data')[key_genes,])))
gene_df$group <- neu$group
gene_df$sample <- neu$sample
write.csv(gene_df, file.path(OUT,"anoikis_human_gene_expr.csv"), row.names=FALSE)

cat("\nDONE 人源 anoikis ->", OUT, "\n")
