#!/usr/bin/env Rscript
# ============================================================
# MEF2C 验证 —— 用质量好的公共数据集（有生物学重复）
# 自有 snRNA n=1/组不可定量，改用 GSE226207(小鼠,15样本) + GSE209552(人源)
# 验证两条链:
#   链A: MEF2C 是否随 TBI 损伤上调（= TGFβRI 下游的间接证据）
#   链B: MEF2C 是否与突触黏附基因(NRXN/CADM/CNTN/CDH) 负相关（= MEF2C 抑制黏附程序）
# 规则: 真实数据, 无模拟; readRDS 后纯度门禁; 伪bulk(按样本)为主、单细胞相关为辅
# ============================================================
suppressMessages({library(Seurat); library(Matrix)})
P <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project"
OUT <- file.path(P, "JNI/analysis/mef2c_verify_out"); dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

gene_ci <- function(mat, g) { h <- rownames(mat)[toupper(rownames(mat))==toupper(g)]; if(length(h)) h[1] else NA }

verify <- function(tag, rds, sample_col, group_cols, adh_genes, astro_mk) {
  cat("\n################ ", tag, " ################\n")
  obj <- readRDS(rds)
  obj <- JoinLayers(obj, assay="RNA")
  cnt <- GetAssayData(obj, assay="RNA", layer="counts")
  dat <- GetAssayData(obj, assay="RNA", layer="data")
  cat("cells:", ncol(obj), " genes:", nrow(obj), "\n")
  # --- 元数据列防御 ---
  md <- colnames(obj@meta.data)
  if (!(sample_col %in% md)) {
    alt <- intersect(c("sample_id","sample","orig.ident","gsm","donor","patient"), md)
    cat("!! sample_col '", sample_col, "' 不存在; 可用候选:", paste(alt,collapse=","), "\n")
    if (length(alt)) sample_col <- alt[1] else { cat("无样本列, 跳过\n"); return(invisible()) }
    cat("   改用 sample_col =", sample_col, "\n")
  }
  group_cols <- intersect(group_cols, md)
  if (length(group_cols)==0) {
    cat("!! 指定的分组列都不存在; 全部 metadata 列:\n"); print(md)
    group_cols <- intersect(c("condition","group","stim","disease","Group","status"), md)
    cat("   改用 group_cols =", paste(group_cols,collapse=","), "\n")
  }

  # --- 纯度门禁 ---
  aqp4 <- gene_ci(dat, astro_mk[1]); snap <- gene_ci(dat, astro_mk[3])
  r <- mean(dat[aqp4,]) / max(mean(dat[snap,]), 1e-6)
  cat(sprintf("纯度: mean(%s)=%.3f  mean(%s)=%.3f  ratio=%.1f\n",
              aqp4, mean(dat[aqp4,]), snap, mean(dat[snap,]), r))
  if (r < 3) cat("  !! 纯度比值<3, 该对象可能不是纯星形胶质——结果谨慎解读\n")

  mef <- gene_ci(dat, "Mef2c")
  if (is.na(mef)) { cat("MEF2C 不在矩阵, 跳过\n"); return(invisible()) }
  adh <- setNames(sapply(adh_genes, function(g) gene_ci(dat,g)), adh_genes)
  adh <- adh[!is.na(adh)]
  cat("MEF2C =", mef, " | 黏附基因命中:", paste(names(adh),collapse=","), "\n")

  # --- 伪bulk: 按样本聚合 counts -> log-CPM ---
  smp <- as.character(obj@meta.data[[sample_col]])
  usmp <- sort(unique(smp))
  pb <- sapply(usmp, function(s) Matrix::rowSums(cnt[, smp==s, drop=FALSE]))
  cpm <- t(t(pb)/colSums(pb))*1e6
  logcpm <- log1p(cpm)
  # 每样本分组信息
  meta_s <- t(sapply(usmp, function(s){
    sub <- obj@meta.data[smp==s,,drop=FALSE]
    sapply(group_cols, function(gc) as.character(sub[[gc]][1]))
  }))
  colnames(meta_s) <- group_cols

  cat("\n--- 链A: MEF2C 伪bulk 表达 (log-CPM) 按样本 ---\n")
  tabA <- data.frame(sample=usmp, meta_s, MEF2C=round(logcpm[mef,usmp],3),
                      n_cell=as.integer(table(smp)[usmp]))
  print(tabA, row.names=FALSE)
  write.csv(tabA, file.path(OUT, paste0(tag,"_MEF2C_pseudobulk.csv")), row.names=FALSE)

  # 分组比较
  for (gc in group_cols) {
    g <- meta_s[,gc]; if (length(unique(g))<2) next
    cat("\n  [", gc, "] MEF2C log-CPM by group:\n")
    for (lv in unique(g)) cat(sprintf("    %-22s n=%d  mean=%.3f\n",
        lv, sum(g==lv), mean(logcpm[mef, usmp[g==lv]])))
    if (length(unique(g))==2) {
      lv <- unique(g)
      tt <- t.test(logcpm[mef,usmp[g==lv[1]]], logcpm[mef,usmp[g==lv[2]]])
      cat(sprintf("    t-test %s vs %s: p=%.4f\n", lv[1], lv[2], tt$p.value))
    }
  }

  cat("\n--- 链B: MEF2C vs 黏附基因 跨样本相关 (Spearman, 伪bulk) ---\n")
  corB <- data.frame()
  for (a in names(adh)) {
    ct <- suppressWarnings(cor.test(logcpm[mef,usmp], logcpm[adh[a],usmp], method="spearman"))
    corB <- rbind(corB, data.frame(adhesion_gene=a, rho=round(ct$estimate,3), p=round(ct$p.value,4)))
  }
  print(corB, row.names=FALSE)
  write.csv(corB, file.path(OUT, paste0(tag,"_MEF2C_adhesion_corr.csv")), row.names=FALSE)

  # --- 单细胞层面: MEF2C-high vs MEF2C-low 星形胶质的黏附基因 ---
  cat("\n--- 链B辅证: 单细胞 MEF2C+ vs MEF2C- 星形胶质 黏附基因均值 ---\n")
  mpos <- dat[mef,] > 0
  cat(sprintf("  MEF2C+ 细胞: %d (%.1f%%)   MEF2C- 细胞: %d\n",
              sum(mpos), 100*mean(mpos), sum(!mpos)))
  scB <- data.frame()
  for (a in names(adh)) {
    hi <- mean(dat[adh[a], mpos]); lo <- mean(dat[adh[a], !mpos])
    wp <- suppressWarnings(wilcox.test(dat[adh[a],mpos], dat[adh[a],!mpos])$p.value)
    scB <- rbind(scB, data.frame(adhesion_gene=a, MEF2Cpos_mean=round(hi,4),
                 MEF2Cneg_mean=round(lo,4), log2FC_pos_vs_neg=round(log2((hi+1e-6)/(lo+1e-6)),3),
                 wilcox_p=signif(wp,3)))
  }
  print(scB, row.names=FALSE)
  write.csv(scB, file.path(OUT, paste0(tag,"_MEF2C_singlecell.csv")), row.names=FALSE)
  invisible()
}

# 黏附基因: 框架塌陷组 NRXN 领衔(NRXN/CNTN/CADM/CDH); NCAM 作阴性对照
adh_mouse <- c("Nrxn1","Nrxn3","Cadm1","Cadm2","Cntn1","Cdh2","Negr1","Ncam1")
adh_human <- c("NRXN1","NRXN3","CADM1","CADM2","CNTN1","CDH2","NEGR1","NCAM1")

verify("GSE226207_mouse",
       file.path(P,"processed_data/human_validation/gse226207_astro_validation.rds"),
       sample_col="sample_id", group_cols=c("timepoint","treatment_group"),
       adh_genes=adh_mouse, astro_mk=c("Aqp4","Gfap","Snap25"))

verify("GSE209552_human",
       file.path(P,"processed_data/human_validation/gse209552_astro_validation.rds"),
       sample_col="sample_id", group_cols=c("timepoint","treatment_group"),
       adh_genes=adh_human, astro_mk=c("AQP4","GFAP","SNAP25"))

cat("\n\nDONE. 输出: ", OUT, "\n")
