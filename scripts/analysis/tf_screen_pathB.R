#!/usr/bin/env Rscript
# ============================================================
# 路 B —— 无偏 TF 筛选: 找"TGFβRI 重编程黏附"的中介转录因子
# 在 3 个可靠公共数据集(有重复)里独立跑, 取跨物种一致的候选
#   GSE226207  小鼠 stab-wound astrocyte (intact vs 3/5dpi)
#   GSE101901+180862 整合 小鼠 astrocyte (control vs tbi)
#   GSE209552  人源 TBI astrocyte (control vs tbi)
# 每个 TF 算两指标:
#   lfc_injury        : 伪bulk 损伤 vs 对照 log2FC  (>0 = 损伤上调)
#   rho_within_injury : 损伤态星形胶质内, TF 与黏附模块的 Spearman 相关
# 候选"塌陷驱动者(抑制型)": lfc>0 且 rho<0  且跨数据集一致
# 候选"丢失的正调控者"    : lfc<0 且 rho>0  且跨数据集一致
# 规则: 真实数据无模拟; readRDS 后纯度门禁; 先全量后推荐
# ============================================================
suppressMessages({library(Seurat); library(Matrix)})
P <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project"
OUT <- file.path(P,"JNI/analysis/tf_screen_out"); dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
tf_mm <- readLines(file.path(P,"scenic_input/databases/allTFs_mm.txt"))
cat("TF 列表:", length(tf_mm), "个(小鼠symbol)\n")
ADH <- c("Nrxn1","Nrxn3","Cadm1","Cadm2","Cntn1","Cdh2")   # 框架 NRXN 领衔塌陷组

gci <- function(rn,g) { h<-rn[toupper(rn)==toupper(g)]; if(length(h)) h[1] else NA }

screen_one <- function(tag, rds, species, astro_filter, sample_col, injured_levels, group_col) {
  cat("\n############", tag, "############\n")
  obj <- readRDS(rds); obj <- JoinLayers(obj, assay="RNA")
  if (!is.null(astro_filter)) {
    keep <- obj@meta.data[[astro_filter$col]] == astro_filter$val
    obj <- obj[, keep]; cat("subset astrocyte:", ncol(obj), "cells\n")
  } else cat("cells(已是astrocyte对象):", ncol(obj), "\n")
  cnt <- GetAssayData(obj,assay="RNA",layer="counts")
  dat <- GetAssayData(obj,assay="RNA",layer="data")
  if (max(dat[,1:min(50,ncol(dat))]) == 0 || identical(dim(dat),dim(cnt)) && max(dat)>50) {
    cat("data层异常, 重新NormalizeData\n"); obj<-NormalizeData(obj,verbose=FALSE)
    dat<-GetAssayData(obj,assay="RNA",layer="data")
  }
  rn <- rownames(dat)
  # 纯度门禁
  aqp4<-gci(rn,"Aqp4"); snap<-gci(rn,"Snap25")
  pr <- mean(dat[aqp4,])/max(mean(dat[snap,]),1e-6)
  cat(sprintf("纯度 Aqp4/Snap25 = %.1f%s\n", pr, ifelse(pr<3," (<3 谨慎)","")))
  # 黏附模块
  ag <- na.omit(sapply(ADH,function(g) gci(rn,g)))
  cat("黏附模块基因:", paste(ag,collapse=","), "\n")
  adh_mod <- Matrix::colMeans(dat[ag,,drop=FALSE])
  # 分组
  grp <- as.character(obj@meta.data[[group_col]])
  inj <- grp %in% injured_levels
  cat("损伤态:", sum(inj), " 对照:", sum(!inj), "\n")
  # TF 集合: 大小写无关匹配(人源对象用小鼠风格symbol) + 损伤态内检出≥5%
  tf_use <- rn[toupper(rn) %in% toupper(tf_mm)]
  cat("TF 命中矩阵:", length(tf_use), "\n")
  detok <- Matrix::rowMeans(cnt[tf_use, inj, drop=FALSE] > 0) >= 0.05
  tf_use <- tf_use[detok]
  cat("筛后可用 TF:", length(tf_use), "(损伤态检出≥5%)\n")
  # 链相关: 损伤态内 TF vs 黏附模块
  Xi <- t(as.matrix(dat[tf_use, inj, drop=FALSE]))
  rho_in <- suppressWarnings(cor(Xi, adh_mod[inj], method="spearman")[,1])
  Xa <- t(as.matrix(dat[tf_use,,drop=FALSE]))
  rho_all <- suppressWarnings(cor(Xa, adh_mod, method="spearman")[,1])
  # 损伤效应: 伪bulk per sample
  smp <- as.character(obj@meta.data[[sample_col]])
  usmp <- unique(smp)
  pb <- sapply(usmp,function(s) Matrix::rowSums(cnt[tf_use,smp==s,drop=FALSE]))
  libsz <- sapply(usmp,function(s) sum(cnt[,smp==s]))
  cpm <- t(t(pb)/libsz)*1e6
  s_inj <- sapply(usmp,function(s) (grp[smp==s][1]) %in% injured_levels)
  lfc <- log2((rowMeans(cpm[,s_inj,drop=FALSE])+1)/(rowMeans(cpm[,!s_inj,drop=FALSE])+1))
  res <- data.frame(TF_upper=toupper(tf_use), TF=tf_use,
                    lfc_injury=round(lfc,3), rho_within=round(rho_in,3),
                    rho_all=round(rho_all,3),
                    pct_detect=round(Matrix::rowMeans(cnt[tf_use,inj,drop=FALSE]>0),3),
                    row.names=NULL)
  write.csv(res, file.path(OUT,paste0(tag,"_TFscreen.csv")), row.names=FALSE)
  cat("输出:", paste0(tag,"_TFscreen.csv"), " (", nrow(res), "TF)\n")
  rm(obj,cnt,dat,Xi,Xa); gc(verbose=FALSE)
  res
}

r1 <- screen_one("GSE226207_ms", file.path(P,"processed_data/human_validation/gse226207_astro_validation.rds"),
  "mouse", NULL, "sample_id", c("3dpi","5dpi"), "timepoint")
r2 <- screen_one("GSE101901_180862_ms", file.path(P,"processed_data/integrated/mouse_integrated_annotated.rds"),
  "mouse", list(col="celltype_auto",val="astrocyte"), "animal", c("tbi"), "condition")
r3 <- screen_one("GSE209552_hs", file.path(P,"processed_data/human_validation/gse209552_astro_validation.rds"),
  "human", NULL, "sample_id", c("tbi"), "treatment_group")

# ---- 跨数据集合并 ----
cat("\n\n################ 跨 3 数据集整合 ################\n")
mg <- Reduce(function(a,b) merge(a,b,by="TF_upper",all=FALSE),
  list(r1[,c("TF_upper","lfc_injury","rho_within")],
       r2[,c("TF_upper","lfc_injury","rho_within")],
       r3[,c("TF_upper","lfc_injury","rho_within")]))
colnames(mg) <- c("TF","lfc_226207","rho_226207","lfc_integ","rho_integ","lfc_209552","rho_209552")
cat("3 数据集共有 TF:", nrow(mg), "\n")
# 一致性打分
lfc <- mg[,c("lfc_226207","lfc_integ","lfc_209552")]
rho <- mg[,c("rho_226207","rho_integ","rho_209552")]
mg$n_repressor   <- rowSums(lfc>0 & rho<0)        # 损伤上调 + 抑制黏附
mg$n_lostact     <- rowSums(lfc<0 & rho>0)        # 损伤下调 + 正调黏附
mg$mean_lfc <- round(rowMeans(lfc),3)
mg$mean_rho <- round(rowMeans(rho),3)
write.csv(mg, file.path(OUT,"CROSS_3datasets_TFscreen.csv"), row.names=FALSE)

cat("\n===== 候选①「塌陷驱动者(损伤上调+抑制黏附)」3 数据集全一致 =====\n")
c1 <- mg[mg$n_repressor==3,]; c1 <- c1[order(c1$mean_rho),]
if(nrow(c1)) print(c1[,c("TF","lfc_226207","lfc_integ","lfc_209552","rho_226207","rho_integ","rho_209552")],row.names=FALSE) else cat("(无 3/3 一致)\n")
cat("\n--- 退一步: 2/3 数据集一致的塌陷驱动者 ---\n")
c1b <- mg[mg$n_repressor==2,]; c1b<-c1b[order(c1b$mean_rho),]
if(nrow(c1b)) print(utils::head(c1b[,c("TF","mean_lfc","mean_rho","n_repressor")],25),row.names=FALSE) else cat("(无)\n")

cat("\n===== 候选②「丢失的正调控者(损伤下调+正调黏附)」3 数据集全一致 =====\n")
c2 <- mg[mg$n_lostact==3,]; c2<-c2[order(-c2$mean_rho),]
if(nrow(c2)) print(utils::head(c2[,c("TF","lfc_226207","lfc_integ","lfc_209552","rho_226207","rho_integ","rho_209552")],25),row.names=FALSE) else cat("(无 3/3 一致)\n")

cat("\nDONE. 输出目录:", OUT, "\n")
