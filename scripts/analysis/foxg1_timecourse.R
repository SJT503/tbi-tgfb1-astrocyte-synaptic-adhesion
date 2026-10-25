#!/usr/bin/env Rscript
# 批注#9 —— FOXG1 在 TBI 后皮层星形胶质细胞的时序表达
# 数据集 GSE271769：重度 CCI 皮层 scRNA，NoTBI/D1/D7/D28 各 n=4（10x h5）
suppressMessages({library(Seurat); library(Matrix)})
set.seed(1)

RAW <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/raw_data/GSE271769"
OUT <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/foxg1_timecourse"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

h5s <- list.files(RAW, pattern="filtered_feature_bc_matrix.h5$", full.names=TRUE)
cat("h5 文件:", length(h5s), "\n")
stopifnot(length(h5s)==16)

meta_tp <- function(f){
  b <- basename(f)
  if (grepl("NoTBI",b)) "NoTBI" else if (grepl("_D1_",b)) "D1" else
  if (grepl("_D7_",b)) "D7" else if (grepl("_D28_",b)) "D28" else NA
}
sid <- function(f) sub("_filtered.*","",basename(f))

# 标志物
M <- list(
  astro   = c("Aqp4","Gfap","Slc1a3","Slc1a2","Aldoc","Gja1"),
  neuron  = c("Snap25","Syt1","Rbfox3","Meg3","Syn1"),
  micro   = c("Cx3cr1","Hexb","C1qa","Ctss","Csf1r"),
  oligo   = c("Plp1","Mbp","Mog","Mobp"),
  opc     = c("Pdgfra","Cspg4"),
  vasc    = c("Cldn5","Flt1","Pdgfrb","Vtn"))

objs <- list()
for (f in h5s){
  m <- Read10X_h5(f)
  if (is.list(m)) m <- m[["Gene Expression"]]
  o <- CreateSeuratObject(m, project=sid(f), min.cells=3, min.features=200)
  o$timepoint <- meta_tp(f); o$sample <- sid(f)
  o[["percent.mt"]] <- PercentageFeatureSet(o, pattern="^mt-")
  o <- subset(o, subset = nFeature_RNA>200 & nFeature_RNA<7000 & percent.mt<20)
  o <- NormalizeData(o, verbose=FALSE)
  objs[[sid(f)]] <- o
  cat(sprintf("  %-26s tp=%-6s cells=%d\n", sid(f), meta_tp(f), ncol(o)))
}

# 合并
obj <- merge(objs[[1]], objs[-1])
obj <- JoinLayers(obj)
cat("合并后总细胞:", ncol(obj), "\n")

# 每细胞按标志物模块打分 → argmax 分类
sc <- sapply(names(M), function(k){
  g <- intersect(M[[k]], rownames(obj))
  Matrix::colMeans(GetAssayData(obj, layer="data")[g,,drop=FALSE])
})
obj$celltype <- colnames(sc)[max.col(sc, ties.method="first")]
cat("\n细胞类型分布:\n"); print(table(obj$celltype))

# 取星形胶质 + 纯度门禁（铁律）
ast <- subset(obj, subset = celltype=="astro")
expr <- GetAssayData(ast, layer="data")
pr <- mean(expr["Aqp4",]) / max(mean(expr["Snap25",]), 1e-6)
cat(sprintf("\n星形胶质子集: n=%d  纯度 Aqp4/Snap25=%.1f\n", ncol(ast), pr))
if (pr < 3) stop("PURITY FAIL: Aqp4 未远高于 Snap25")

# Foxg1 在星形胶质细胞中：每样本伪bulk均值（lognorm）
stopifnot("Foxg1" %in% rownames(ast))
fox <- GetAssayData(ast, layer="data")["Foxg1",]
df <- data.frame(sample=ast$sample, timepoint=ast$timepoint, foxg1=fox)
agg <- aggregate(foxg1~sample+timepoint, df, mean)
agg$pct <- aggregate(foxg1~sample+timepoint, df, function(x) mean(x>0))$foxg1
agg$timepoint <- factor(agg$timepoint, levels=c("NoTBI","D1","D7","D28"))
agg <- agg[order(agg$timepoint),]
write.csv(agg, file.path(OUT,"foxg1_astro_per_sample.csv"), row.names=FALSE)

cat("\n=== FOXG1 星形胶质细胞表达（每样本伪bulk均值）===\n")
print(agg)

# 各时间点汇总 + 统计（每组 n=4）
summ <- aggregate(foxg1~timepoint, agg, function(x) c(mean=mean(x), se=sd(x)/sqrt(length(x)), n=length(x)))
summ <- do.call(data.frame, summ)
cat("\n=== 各时间点汇总 ===\n"); print(summ)

kw <- kruskal.test(foxg1~timepoint, agg)
cat(sprintf("\nKruskal-Wallis（时间点间差异）: p=%.4f\n", kw$p.value))
base <- agg$foxg1[agg$timepoint=="NoTBI"]
for (tp in c("D1","D7","D28")){
  v <- agg$foxg1[agg$timepoint==tp]
  wt <- wilcox.test(v, base)
  cat(sprintf("  %s vs NoTBI: 均值 %.3f vs %.3f  (%.0f%%)  Wilcoxon p=%.3f\n",
      tp, mean(v), mean(base), 100*mean(v)/mean(base), wt$p.value))
}
write.csv(summ, file.path(OUT,"foxg1_timecourse_summary.csv"), row.names=FALSE)
saveRDS(ast, file.path(OUT,"GSE271769_astro.rds"))
cat("\nDONE ->", OUT, "\n")
