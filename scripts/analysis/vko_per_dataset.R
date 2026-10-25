#!/usr/bin/env Rscript
# 批注#5/#17 —— 逐数据集 FOXG1 虚拟敲除 + GSE101901/GSE180862 拆分
# 在每个独立公共数据集的星形胶质细胞上：
#   (1) FOXG1 损伤效应 lfc + 与黏附模块相关 rho（拆开统计，不再合并）
#   (2) Ridge-GRN 虚拟敲除 FOXG1 → 黏附程序基因预测位移（CellOracle 核心算法）
# 规则：真实数据；readRDS 后纯度门禁；先全量后汇总
suppressMessages({library(Seurat); library(Matrix)})
set.seed(1)
P  <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project"
OUT<- file.path(P,"JNI/analysis/vko_per_dataset"); dir.create(OUT,showWarnings=FALSE,recursive=TRUE)
tf_mm <- toupper(readLines(file.path(P,"scenic_input/databases/allTFs_mm.txt")))
ADH   <- c("NRXN1","NRXN3","CADM1","CADM2","CNTN1","CDH2")

gci <- function(rn,g){ h<-rn[toupper(rn)==toupper(g)]; if(length(h)) h[1] else NA }

# ---- 单数据集分析 ----
run_ds <- function(tag, dat, grp, injured){
  rn <- rownames(dat)
  aqp4<-gci(rn,"Aqp4"); snap<-gci(rn,"Snap25"); foxg<-gci(rn,"Foxg1")
  pr <- mean(dat[aqp4,])/max(mean(dat[snap,]),1e-6)
  cat(sprintf("\n#### %s  星形胶质 n=%d  纯度Aqp4/Snap25=%.1f\n", tag, ncol(dat), pr))
  if(is.na(foxg)){ cat("  FOXG1 不在矩阵，跳过\n"); return(NULL) }
  inj <- grp %in% injured
  cat(sprintf("  损伤态 %d  对照 %d\n", sum(inj), sum(!inj)))
  adg <- na.omit(sapply(ADH, function(g) gci(rn,g)))
  adh_mod <- Matrix::colMeans(dat[adg,,drop=FALSE])
  # (1) FOXG1 损伤效应 + 链相关
  lfc <- log2((mean(dat[foxg,inj])+0.01)/(mean(dat[foxg,!inj])+0.01))
  rho <- suppressWarnings(cor(dat[foxg,inj], adh_mod[inj], method="spearman"))
  # (2) Ridge-GRN 虚拟敲除
  det <- Matrix::rowMeans(dat>0) >= 0.05
  g_use <- rn[det]
  cells <- 1:ncol(dat)
  if(length(cells)>8000) cells <- sample(cells,8000)
  X <- t(as.matrix(dat[g_use,cells,drop=FALSE]))      # cell x gene
  mu<-colMeans(X); sd<-apply(X,2,sd); sd[sd==0]<-1
  Z <- scale(X,center=mu,scale=sd)
  tfs <- intersect(g_use, rn[toupper(rn)%in%tf_mm])
  tv  <- apply(Z[,tfs,drop=FALSE],2,var)
  tfs <- names(sort(tv,decreasing=TRUE))[1:min(150,length(tfs))]
  if(!(foxg%in%tfs)) tfs <- c(tfs[-length(tfs)],foxg)
  Xtf <- Z[,tfs,drop=FALSE]
  B <- solve(crossprod(Xtf)+diag(length(tfs)), crossprod(Xtf,Z))   # nTF x nGene
  fp <- which(tfs==foxg)
  dF <- ((0-mu[foxg])/sd[foxg]) - Z[,foxg]
  s1 <- outer(dF, B[fp,])
  s2 <- s1 + s1[,tfs,drop=FALSE] %*% B
  ms <- colMeans(s2); names(ms)<-colnames(Z)
  ash <- ms[adg]
  u <- suppressWarnings(wilcox.test(ash, ms, alternative="less"))
  cat(sprintf("  FOXG1: lfc_injury=%+.3f  rho_adh=%+.3f\n", lfc, rho))
  cat("  虚拟敲除黏附基因位移:\n")
  for(g in adg) cat(sprintf("    %-8s %+.4f\n", g, ms[g]))
  cat(sprintf("  黏附均值=%+.4f  背景中位=%+.4f  下调 %d/%d  Wilcoxon(黏附更低) p=%.2e\n",
      mean(ash), median(ms), sum(ash<0), length(adg), u$p.value))
  data.frame(dataset=tag, n_astro=ncol(dat), purity=round(pr,1),
             foxg1_lfc_injury=round(lfc,3), foxg1_rho_adh=round(rho,3),
             ko_adh_mean_shift=round(mean(ash),4), ko_adh_down=sprintf("%d/%d",sum(ash<0),length(adg)),
             ko_wilcox_p=signif(u$p.value,3))
}

getdat <- function(o){
  o <- JoinLayers(o, assay="RNA")
  d <- GetAssayData(o, assay="RNA", layer="data")
  if(max(d[,1:min(30,ncol(d))])==0){ o<-NormalizeData(o,verbose=FALSE); d<-GetAssayData(o,assay="RNA",layer="data") }
  d
}

res <- list()

## --- GSE101901 / GSE180862：来自整合对象，按 dataset 列拆 ---
cat("载入整合对象 ...\n")
integ <- readRDS(file.path(P,"processed_data/integrated/mouse_integrated_annotated.rds"))
ast_i <- integ[, integ$celltype_auto=="astrocyte"]
for(ds in c("GSE101901","GSE180862")){
  sub <- ast_i[, ast_i$dataset==ds]
  d <- getdat(sub)
  res[[ds]] <- run_ds(ds, d, as.character(sub$condition), "tbi")
}
rm(integ,ast_i); gc(verbose=FALSE)

## --- GSE226207（小鼠 stab-wound）---
o <- readRDS(file.path(P,"processed_data/human_validation/gse226207_astro_validation.rds"))
res[["GSE226207"]] <- run_ds("GSE226207", getdat(o), as.character(o$timepoint), c("3dpi","5dpi"))

## --- GSE209552（人源 TBI）---
o <- readRDS(file.path(P,"processed_data/human_validation/gse209552_astro_validation.rds"))
res[["GSE209552"]] <- run_ds("GSE209552", getdat(o), as.character(o$treatment_group), "tbi")

## --- GSE271769（重度CCI皮层，时序）---
o <- readRDS(file.path(P,"JNI/analysis/foxg1_timecourse/GSE271769_astro.rds"))
res[["GSE271769"]] <- run_ds("GSE271769", getdat(o), as.character(o$timepoint), c("D1","D7","D28"))

tab <- do.call(rbind, res)
write.csv(tab, file.path(OUT,"vko_per_dataset_summary.csv"), row.names=FALSE)
cat("\n\n================ 逐数据集汇总 ================\n")
print(tab, row.names=FALSE)
cat("\nDONE ->", OUT, "\n")
