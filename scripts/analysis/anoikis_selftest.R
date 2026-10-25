#!/usr/bin/env Rscript
# 自测 snRNA-seq 神经元 anoikis signature (3 组: Sham / CCI_Scr / CCI_TGFbetaR1_KO)
# 注: 三组各为单次测序, 定量结果有 batch 局限, 但 3 组方向一致性可作为机制佐证
suppressMessages({library(Seurat); library(Matrix)})
set.seed(1)
P <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project"
OUT <- file.path(P,"JNI/analysis/anoikis_out")

# 合并 excitatory + inhibitory 神经元
exc <- readRDS(file.path(P,"wet_lab_data/results/seurat/excitatory_v2_unified.rds"))
inh <- readRDS(file.path(P,"wet_lab_data/results/seurat/inhibitory_v2_unified.rds"))
DefaultAssay(exc)<-"RNA"; DefaultAssay(inh)<-"RNA"
# 兼容 Seurat v4 (Assay class) 与 v5 (Assay5)
safe_join <- function(o){
  if (inherits(o[["RNA"]], "Assay5")) tryCatch(JoinLayers(o, assay="RNA"), error=function(e) o)
  else o
}
exc <- safe_join(exc); inh <- safe_join(inh)
neu <- merge(exc, inh)
neu <- safe_join(neu)
cat("神经元 (Exc+Inh):", ncol(neu), "\n")

# 纯度门禁
expr <- GetAssayData(neu, layer="data")
pr <- mean(expr["Snap25",])/max(mean(expr["Aqp4",]),1e-6)
cat(sprintf("纯度 Snap25/Aqp4=%.1f\n", pr))
if (pr < 3) stop("PURITY FAIL")

# 标准化(防止 layer 异常)
if (max(GetAssayData(neu, layer="data")[,1:min(50,ncol(neu))]) == 0)
  neu <- NormalizeData(neu, verbose=FALSE)

survival_genes <- c("Ptk2","Akt1","Akt2","Akt3","Bcl2","Mcl1","Bcl2l1","Birc5","Mapk1")
death_genes    <- c("Bcl2l11","Bmf","Bax","Bad","Casp9","Casp3","Cycs","Apaf1")

neu <- AddModuleScore(neu, features=list(survival_genes), name="survival_score", search=TRUE)
neu <- AddModuleScore(neu, features=list(death_genes),    name="death_score",    search=TRUE)
neu$survival <- neu$survival_score1
neu$death    <- neu$death_score1
neu$anoikis  <- neu$death - neu$survival

cat("\n分组:\n"); print(table(neu$orig.ident))

cell_df <- data.frame(sample=neu$orig.ident, group=neu$orig.ident,
                      survival=neu$survival, death=neu$death, anoikis=neu$anoikis)
# 标准化分组名
cell_df$group <- factor(cell_df$group, levels=c("Sham","CCI_Scr","CCI_TGFbetaR1_KO"),
                        labels=c("Sham","CCI","CCI+TGFβRI-KO"))
write.csv(cell_df, file.path(OUT,"anoikis_selftest_per_cell.csv"), row.names=FALSE)

cat("\n=== 各组 anoikis 评分 (每细胞均值±SE) ===\n")
for (s in c("survival","death","anoikis")){
  m <- tapply(cell_df[[s]], cell_df$group, mean)
  se<- tapply(cell_df[[s]], cell_df$group, function(x) sd(x)/sqrt(length(x)))
  cat(s,":\n"); for (g in names(m)) cat(sprintf("  %-22s %+.3f ± %.3f\n",g,m[g],se[g]))
  p_sham_cci <- wilcox.test(cell_df[[s]][cell_df$group=="Sham"], cell_df[[s]][cell_df$group=="CCI"])$p.value
  p_cci_ko   <- wilcox.test(cell_df[[s]][cell_df$group=="CCI"], cell_df[[s]][cell_df$group=="CCI+TGFβRI-KO"])$p.value
  cat(sprintf("    Sham vs CCI  p=%.2e\n", p_sham_cci))
  cat(sprintf("    CCI  vs KO   p=%.2e\n", p_cci_ko))
}

# 基因表达 for heatmap
key_genes <- intersect(c(death_genes, survival_genes), rownames(neu))
gene_df <- data.frame(t(as.matrix(GetAssayData(neu, layer='data')[key_genes,])))
gene_df$group <- as.character(cell_df$group)
write.csv(gene_df, file.path(OUT,"anoikis_selftest_gene_expr.csv"), row.names=FALSE)
cat("DONE 自测 anoikis ->", OUT, "\n")
