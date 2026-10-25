#!/usr/bin/env Rscript
# 人源 GSE209552 raw mtx 处理 → 神经元 anoikis 分析
suppressMessages({library(Seurat); library(Matrix)})
set.seed(1)
RAW <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/raw_data/GSE209552/extracted_all"
OUT <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/anoikis_out"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

# 找所有 mtx 三件套样本
mtx_files <- list.files(RAW, pattern="_matrix\\.mtx\\.gz$", full.names=TRUE)
cat("发现 mtx 样本:", length(mtx_files), "\n")

# 分组: HuBrainCTL 开头 = CTL; MJ_TBI 开头 = TBI
get_group <- function(f){
  if (grepl("HuBrainCTL", f)) "CTL" else if (grepl("MJ_TBI", f)) "TBI" else NA
}
get_sample <- function(f) sub("_matrix\\.mtx\\.gz$","",basename(f))

objs <- list()
for (mf in mtx_files){
  prefix <- sub("_matrix\\.mtx\\.gz$","",mf)
  bf <- paste0(prefix, "_barcodes.tsv.gz")
  ff <- paste0(prefix, "_features.tsv.gz")
  if (!file.exists(bf) || !file.exists(ff)) { cat("  缺三件套:", basename(prefix), "\n"); next }
  grp <- get_group(mf)
  if (is.na(grp)) { cat("  分组未知:", basename(mf), "\n"); next }
  sid <- get_sample(mf)
  m <- readMM(mf)
  bc <- readLines(bf)
  feat <- read.table(ff, sep="\t", header=FALSE)
  if (ncol(feat) >= 2) rownames(m) <- make.unique(feat[,2]) else rownames(m) <- make.unique(feat[,1])
  colnames(m) <- paste0(sid, "_", bc)
  cat(sprintf("  %s  %s  %dgene x %dcell\n", grp, sid, nrow(m), ncol(m)))
  o <- CreateSeuratObject(m, project=sid, min.cells=3, min.features=200)
  if (ncol(o) < 50) next
  o$group <- grp; o$sample <- sid
  o[["percent.mt"]] <- PercentageFeatureSet(o, pattern="^MT-")
  o <- subset(o, subset = nFeature_RNA>200 & percent.mt<25)
  o <- NormalizeData(o, verbose=FALSE)
  objs[[sid]] <- o
}
cat("\n有效样本:", length(objs), "\n")
obj <- merge(objs[[1]], objs[-1])
obj <- JoinLayers(obj, assay="RNA")
cat("合并细胞:", ncol(obj), "  分组:\n"); print(table(obj$group))

# 标志物 (人源大写)
M <- list(
  neuron  = c("SNAP25","SYT1","RBFOX3","MEG3","SYN1","STMN2"),
  astro   = c("AQP4","GFAP","SLC1A3","SLC1A2","ALDH1L1","GJA1"),
  micro   = c("CX3CR1","HEXB","C1QA","CSF1R","P2RY12"),
  oligo   = c("PLP1","MBP","MOG","MOBP"),
  opc     = c("PDGFRA","CSPG4"),
  vasc    = c("CLDN5","FLT1","PDGFRB"))

sc <- sapply(names(M), function(k){
  g <- intersect(M[[k]], rownames(obj))
  if (length(g)==0) return(rep(0, ncol(obj)))
  Matrix::colMeans(GetAssayData(obj, layer="data")[g,,drop=FALSE])
})
obj$celltype <- colnames(sc)[max.col(sc, ties.method="first")]
cat("\n细胞类型:\n"); print(table(obj$celltype))

# 神经元 + 纯度
neu <- subset(obj, subset=celltype=="neuron")
expr <- GetAssayData(neu, layer="data")
pr <- mean(expr["SNAP25",])/max(mean(expr["AQP4",]),1e-6)
cat(sprintf("\n神经元: n=%d  纯度 SNAP25/AQP4=%.1f\n", ncol(neu), pr))
if (pr < 3) stop("PURITY FAIL")

survival_genes <- c("PTK2","AKT1","AKT2","AKT3","BCL2","MCL1","BCL2L1","BIRC5","MAPK1")
death_genes    <- c("BCL2L11","BMF","BAX","BAD","CASP9","CASP3","CYCS","APAF1")

neu <- AddModuleScore(neu, features=list(survival_genes), name="survival_score", search=TRUE)
neu <- AddModuleScore(neu, features=list(death_genes),    name="death_score",    search=TRUE)
neu$survival <- neu$survival_score1
neu$death    <- neu$death_score1
neu$anoikis  <- neu$death - neu$survival

cell_df <- data.frame(sample=neu$sample, group=neu$group,
                      survival=neu$survival, death=neu$death, anoikis=neu$anoikis)
write.csv(cell_df, file.path(OUT,"anoikis_human_per_cell.csv"), row.names=FALSE)

# 样本伪bulk
agg <- aggregate(cbind(survival,death,anoikis) ~ sample+group, cell_df, mean)
write.csv(agg, file.path(OUT,"anoikis_human_per_sample.csv"), row.names=FALSE)

cat("\n=== 人源 anoikis 评分 (每组均值±SE) ===\n")
for (s in c("survival","death","anoikis")){
  m <- tapply(agg[[s]], agg$group, mean)
  se<- tapply(agg[[s]], agg$group, function(x) sd(x)/sqrt(length(x)))
  for (g in names(m)) cat(sprintf("  %s %-4s %+.3f±%.3f (n=%d)\n", s, g, m[g], se[g], sum(agg$group==g)))
  v1<-agg[[s]][agg$group=="CTL"]; v2<-agg[[s]][agg$group=="TBI"]
  p <- wilcox.test(v1, v2)$p.value
  cat(sprintf("    CTL vs TBI: p=%.3e\n", p))
}

# 基因表达 for heatmap
key_genes <- intersect(c(death_genes, survival_genes), rownames(neu))
gene_df <- data.frame(t(as.matrix(GetAssayData(neu, layer='data')[key_genes,])))
gene_df$group <- neu$group
gene_df$sample <- neu$sample
write.csv(gene_df, file.path(OUT,"anoikis_human_gene_expr.csv"), row.names=FALSE)
cat("DONE 人源 ->", OUT, "\n")
