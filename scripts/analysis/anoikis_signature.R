#!/usr/bin/env Rscript
# 图7 panel C: 神经元 anoikis signature 时序变化 (GSE271769 重度CCI皮层数据)
# 验证: TBI 后皮层神经元呈现"失巢凋亡样"转录指纹——
#   存活通路 (FAK/Akt) DOWN, 凋亡执行通路 (caspase/Bim/Bmf) UP
suppressMessages({library(Seurat); library(Matrix)})
set.seed(1)

RAW <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/raw_data/GSE271769"
OUT <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/anoikis_out"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

h5s <- list.files(RAW, pattern="filtered_feature_bc_matrix.h5$", full.names=TRUE)
stopifnot(length(h5s)==16)

meta_tp <- function(f){
  b <- basename(f)
  if (grepl("NoTBI",b)) "NoTBI" else if (grepl("_D1_",b)) "D1" else
  if (grepl("_D7_",b)) "D7" else if (grepl("_D28_",b)) "D28" else NA
}
sid <- function(f) sub("_filtered.*","",basename(f))

# 标志物
M <- list(
  neuron  = c("Snap25","Syt1","Rbfox3","Meg3","Syn1"),
  astro   = c("Aqp4","Gfap","Slc1a3","Slc1a2","Aldoc","Gja1"),
  micro   = c("Cx3cr1","Hexb","C1qa","Ctss","Csf1r"),
  oligo   = c("Plp1","Mbp","Mog","Mobp"),
  opc     = c("Pdgfra","Cspg4"),
  vasc    = c("Cldn5","Flt1","Pdgfrb","Vtn"))

# anoikis 基因签名 (文献综合)
survival_genes  <- c("Ptk2","Akt1","Akt2","Akt3","Bcl2","Mcl1","Bcl2l1","Birc5","Mapk1")  # FAK/Akt 存活通路
death_genes     <- c("Bcl2l11","Bmf","Bax","Bad","Casp9","Casp3","Cycs","Apaf1")           # 凋亡执行通路

objs <- list()
for (f in h5s){
  m <- Read10X_h5(f); if (is.list(m)) m <- m[["Gene Expression"]]
  o <- CreateSeuratObject(m, project=sid(f), min.cells=3, min.features=200)
  o$timepoint <- meta_tp(f); o$sample <- sid(f)
  o[["percent.mt"]] <- PercentageFeatureSet(o, pattern="^mt-")
  o <- subset(o, subset = nFeature_RNA>200 & nFeature_RNA<7000 & percent.mt<20)
  o <- NormalizeData(o, verbose=FALSE)
  objs[[sid(f)]] <- o
}
obj <- merge(objs[[1]], objs[-1]); obj <- JoinLayers(obj)
cat("合并总细胞:", ncol(obj), "\n")

# 细胞分类
sc <- sapply(names(M), function(k){
  g <- intersect(M[[k]], rownames(obj))
  Matrix::colMeans(GetAssayData(obj, layer="data")[g,,drop=FALSE])
})
obj$celltype <- colnames(sc)[max.col(sc, ties.method="first")]
cat("细胞类型:\n"); print(table(obj$celltype))

# 神经元子集 + 纯度门禁
neu <- subset(obj, subset=celltype=="neuron")
expr <- GetAssayData(neu, layer="data")
pr <- mean(expr["Snap25",])/max(mean(expr["Aqp4",]),1e-6)
cat(sprintf("\n神经元子集: n=%d  纯度 Snap25/Aqp4=%.1f\n", ncol(neu), pr))
if (pr < 3) stop("PURITY FAIL: neuron Snap25 not >> Aqp4")

# 计算 anoikis signature score
neu <- AddModuleScore(neu, features=list(survival_genes), name="survival_score", search=TRUE)
neu <- AddModuleScore(neu, features=list(death_genes),    name="death_score",    search=TRUE)
neu$survival <- neu$survival_score1
neu$death    <- neu$death_score1
neu$anoikis  <- neu$death - neu$survival   # 越大越倾向 anoikis

# 每样本伪 bulk (平均)
df <- data.frame(sample=neu$sample, timepoint=neu$timepoint,
                 survival=neu$survival, death=neu$death, anoikis=neu$anoikis)
agg <- aggregate(cbind(survival,death,anoikis) ~ sample+timepoint, df, mean)
agg$timepoint <- factor(agg$timepoint, levels=c("NoTBI","D1","D7","D28"))
agg <- agg[order(agg$timepoint),]
write.csv(agg, file.path(OUT,"anoikis_per_sample.csv"), row.names=FALSE)

cat("\n=== 神经元 anoikis 时序 (每样本) ===\n")
print(agg, row.names=FALSE)

cat("\n=== 各时间点汇总 (均值±SE) ===\n")
for (s in c("survival","death","anoikis")){
  m <- tapply(agg[[s]], agg$timepoint, mean)
  se<- tapply(agg[[s]], agg$timepoint, function(x) sd(x)/sqrt(length(x)))
  cat(s,":\n"); for (tp in names(m)) cat(sprintf("  %-6s  %+.3f ± %.3f\n",tp,m[tp],se[tp]))
  # 各 TBI 时间点 vs NoTBI
  base <- agg[[s]][agg$timepoint=="NoTBI"]
  for (tp in c("D1","D7","D28")){
    v <- agg[[s]][agg$timepoint==tp]
    p <- wilcox.test(v, base)$p.value
    cat(sprintf("    %s vs NoTBI: p=%.3f\n", tp, p))
  }
}

# Kruskal 整体
for (s in c("survival","death","anoikis")){
  p <- kruskal.test(agg[[s]] ~ agg$timepoint)$p.value
  cat(sprintf("\nKruskal %s ~ timepoint: p=%.3f\n", s, p))
}
cat("\nDONE ->", OUT, "\n")
