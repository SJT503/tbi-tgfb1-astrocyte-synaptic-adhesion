# Script 104b: 谷氨酸面板 v2 — 剩余 7 数据集（104 在 GSE290150 卡死后的修复版）
# ═══════════════════════════════════════════════════════════════════
# 104 已出 4 个：GSE101901/180862_Cor/180862_Hip/247339（结果见 104 日志）
# 本脚本只跑剩余 7 个，并修两个问题：
#   ① find_astro 加 4 万细胞上限（聚类前下采样）防大对象卡死
#   ② GSE290150 先在稀疏矩阵按细胞行子集，再转置（避免全矩阵操作）
# 每步加进度打印

suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(dplyr) })
options(future.globals.maxSize=4e9)

BASE <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project"
RAW  <- file.path(BASE, "raw_data")
HUM  <- file.path(BASE, "wet_lab_data/results/human_validation_v2/processed.rds")
FIGOUT <- file.path(BASE, "wet_lab_data/results/figures/Fig_Homeostatic_Deepdive")

MK <- list(Neuron=c("Snap25","Syt1","Rbfox3"), Astro=c("Aqp4","Gfap","Slc1a2"),
  Oligo=c("Plp1","Mbp","Mobp"), Micro=c("Csf1r","P2ry12","Cx3cr1"),
  OPC=c("Pdgfra","Cspg4"), Endo=c("Cldn5","Flt1"))

find_astro <- function(seu, tag="") {
  if (ncol(seu) > 40000) {
    set.seed(1); keep <- sample(colnames(seu), 40000)
    seu <- seu[, keep]
    cat(sprintf("    [%s] 细胞>4万，下采样到 40000 再聚类\n", tag))
  }
  cat(sprintf("    [%s] 聚类中 (n=%d)...\n", tag, ncol(seu))); flush.console()
  seu <- NormalizeData(seu, verbose=FALSE)
  seu <- FindVariableFeatures(seu, nfeatures=2000, verbose=FALSE)
  seu <- ScaleData(seu, verbose=FALSE); seu <- RunPCA(seu, npcs=30, verbose=FALSE)
  seu <- FindNeighbors(seu, dims=1:30, verbose=FALSE)
  seu <- FindClusters(seu, resolution=0.5, verbose=FALSE)
  ex <- GetAssayData(seu, layer="data", assay="RNA")
  cls <- sort(unique(seu$seurat_clusters))
  mk <- if (toupper(rownames(ex)[1])==rownames(ex)[1]) lapply(MK,toupper) else MK
  sc <- sapply(names(mk), function(ct) {
    g <- intersect(mk[[ct]], rownames(ex))
    sapply(cls, function(cl) mean(as.matrix(ex[g, colnames(seu)[seu$seurat_clusters==cl], drop=FALSE]))) })
  acl <- cls[colnames(sc)[apply(sc,1,which.max)]=="Astro"]
  seu[, seu$seurat_clusters %in% acl]
}

glu_panel <- function(astro, g_tbi, g_ctrl, dataset, model, severity, region, platform) {
  upper <- toupper(rownames(astro)[1])==rownames(astro)[1]
  GLU <- if (upper) c("GLUL","SLC1A2","SLC1A3") else c("Glul","Slc1a2","Slc1a3")
  grp <- astro$Group
  n1 <- sum(grp==g_tbi,na.rm=TRUE); n0 <- sum(grp==g_ctrl,na.rm=TRUE)
  if (n1<15 || n0<15) { cat(sprintf("  [%s] 星形太少(%d/%d) 跳过\n",dataset,n1,n0)); return(NULL) }
  counts <- GetAssayData(astro, layer="counts", assay="RNA")
  med <- c(median(colSums(counts[,colnames(astro)[grp==g_tbi],drop=FALSE])),
           median(colSums(counts[,colnames(astro)[grp==g_ctrl],drop=FALSE])))
  set.seed(42)
  ds <- SampleUMI(counts, max.umi=round(min(med)*0.9), upsample=FALSE, verbose=FALSE)
  rownames(ds) <- rownames(counts); colnames(ds) <- colnames(counts)
  d <- CreateSeuratObject(counts=ds); d$Group <- as.character(grp)
  d <- NormalizeData(d, verbose=FALSE)
  de <- GetAssayData(d, layer="data", assay="RNA")
  row <- data.frame(dataset=dataset, model=model, severity=severity, region=region,
                    platform=platform, n_astro=ncol(astro), n_TBI=n1, n_ctrl=n0)
  for (g in GLU) {
    gg <- toupper(g)
    if (!g %in% rownames(de)) { row[[paste0(gg,"_pct")]] <- NA; row[[paste0(gg,"_p")]] <- NA; next }
    vt <- as.numeric(de[g, grp==g_tbi]); vc <- as.numeric(de[g, grp==g_ctrl])
    row[[paste0(gg,"_pct")]] <- round((mean(vt)-mean(vc))/(mean(vc)+1e-9)*100,1)
    row[[paste0(gg,"_p")]] <- signif(tryCatch(wilcox.test(vt,vc)$p.value,error=function(e)NA),3)
  }
  cat(sprintf("  [%s] 星形 %d (T%d/C%d) | GLUL %+.1f%% / SLC1A2 %+.1f%% / SLC1A3 %+.1f%%\n",
              dataset, ncol(astro), n1, n0, row$GLUL_pct, row$SLC1A2_pct, row$SLC1A3_pct))
  row
}

read_gsm <- function(d, gsm) {
  f <- list.files(d, pattern=paste0("^",gsm,"[_.]"), full.names=TRUE)
  mtx <- grep("matrix\\.mtx", f, value=TRUE)
  bc  <- grep("barcodes\\.tsv", f, value=TRUE)
  ft  <- grep("(features|genes)\\.tsv", f, value=TRUE)
  if (length(mtx)!=1||length(bc)!=1||length(ft)!=1) stop(paste("file mismatch:",gsm))
  nc <- length(strsplit(readLines(gzfile(ft),n=1),"\t")[[1]])
  ReadMtx(mtx=mtx, cells=bc, features=ft, feature.column=min(nc,2))
}
build_merged <- function(d, gsm_grp, mt="^mt-") {
  objs <- lapply(names(gsm_grp), function(gsm) {
    o <- CreateSeuratObject(read_gsm(d,gsm), min.cells=3, min.features=200)
    o$Group <- gsm_grp[[gsm]]; o })
  seu <- if (length(objs)==1) objs[[1]] else merge(objs[[1]], objs[-1])
  if (length(objs)>1) seu <- JoinLayers(seu)
  seu[["pm"]] <- PercentageFeatureSet(seu, pattern=mt)
  subset(seu, nFeature_RNA>=200 & pm<25)
}

RES <- list()

# ── GSE290150 中度CCI皮层 ──
cat("\n== GSE290150 (中度CCI/皮层/scRNA) ==\n"); flush.console()
tryCatch({ d <- file.path(RAW,"GSE290150")
  cat("  读矩阵(fread 快读 2.66亿三元组)...\n"); flush.console()
  hdr <- scan(gzfile(file.path(d,"GSE290150_matrix.mtx.gz")), what=integer(), skip=1, nmax=3, quiet=TRUE)
  dt <- data.table::fread(file.path(d,"GSE290150_matrix.mtx.gz"), skip=2, header=FALSE)
  bc <- read.table(gzfile(file.path(d,"GSE290150_features.tsv.gz")),sep="\t")[,1]  # 实为 barcodes
  gn <- read.table(gzfile(file.path(d,"GSE290150_barcodes.tsv.gz")),sep="\t")[,1]  # 实为 genes
  cat(sprintf("  三元组 %d 行，构建稀疏矩阵 %d x %d...\n",nrow(dt),hdr[1],hdr[2])); flush.console()
  m <- sparseMatrix(i=dt$V1, j=dt$V2, x=dt$V3, dims=hdr[1:2])  # cell x gene
  rm(dt); gc(verbose=FALSE)
  grpf <- sub("^[^.]*\\.","",bc)
  keep <- grepl("Normothermia",grpf)
  cat(sprintf("  Normothermia 细胞 %d / 总 %d，子集后转置...\n",sum(keep),length(bc))); flush.console()
  m <- t(m[keep,])                       # 按细胞行子集 → 转置 gene x cell
  rownames(m) <- make.unique(as.character(gn)); colnames(m) <- bc[keep]
  seu <- CreateSeuratObject(counts=m, min.cells=3, min.features=200)
  seu$Group <- ifelse(grepl("^Sham",grpf[keep][match(colnames(seu),bc[keep])]),"Sham","TBI")
  seu[["pm"]] <- PercentageFeatureSet(seu,pattern="^mt-")
  seu <- subset(seu, nFeature_RNA>=200 & pm<25)
  RES[["GSE290150"]] <- glu_panel(find_astro(seu,"GSE290150"),"TBI","Sham","GSE290150","CCI","中度","皮层","scRNA")
  rm(m,seu); gc(verbose=FALSE) }, error=function(e) cat("  GSE290150 失败:",conditionMessage(e),"\n"))

# ── GSE226207 皮层刺伤 ──
cat("\n== GSE226207 (皮层刺伤/scRNA) ==\n"); flush.console()
tryCatch({ gg <- as.list(c(
    GSM7068147="Sham",GSM7068152="Sham",GSM7068156="Sham",GSM7068157="Sham",GSM7068159="Sham",
    GSM7068148="TBI", GSM7068149="TBI", GSM7068153="TBI", GSM7068158="TBI", GSM7068160="TBI"))
  seu <- build_merged(file.path(RAW,"GSE226207"), gg)
  RES[["GSE226207"]] <- glu_panel(find_astro(seu,"GSE226207"),"TBI","Sham","GSE226207","刺伤","刺伤","皮层","scRNA")
  rm(seu); gc(verbose=FALSE) }, error=function(e) cat("  GSE226207 失败:",conditionMessage(e),"\n"))

# ── GSE271769 重度CCI皮层 (10x h5) ──
cat("\n== GSE271769 (重度CCI/皮层/scRNA h5) ==\n"); flush.console()
tryCatch({ d <- file.path(RAW,"GSE271769")
  h5 <- list.files(d, pattern="\\.h5$", full.names=TRUE)
  objs <- lapply(h5, function(f) {
    grp <- if (grepl("NoTBI",basename(f))) "Sham" else "TBI"
    o <- CreateSeuratObject(Read10X_h5(f), min.cells=3, min.features=200)
    o$Group <- grp; o })
  seu <- merge(objs[[1]],objs[-1]); seu <- JoinLayers(seu)
  seu[["pm"]] <- PercentageFeatureSet(seu,pattern="^mt-")
  seu <- subset(seu, nFeature_RNA>=200 & pm<25)
  RES[["GSE271769"]] <- glu_panel(find_astro(seu,"GSE271769"),"TBI","Sham","GSE271769","CCI","重度","皮层","scRNA")
  rm(objs,seu); gc(verbose=FALSE) }, error=function(e) cat("  GSE271769 失败:",conditionMessage(e),"\n"))

# ── GSE277487 重度CCI皮层 (snRNA) ──
cat("\n== GSE277487 (重度CCI/皮层/snRNA) ==\n"); flush.console()
tryCatch({ d <- file.path(RAW,"GSE277487")
  gsms <- unique(sub("_.*$","",list.files(d,pattern="matrix\\.mtx")))
  gg <- setNames(lapply(gsms, function(x) {
    fn <- list.files(d,pattern=paste0("^",x,"_"))[1]
    if (grepl("Ctrl",fn)) "Sham" else "TBI" }), gsms)
  seu <- build_merged(d, gg)
  RES[["GSE277487"]] <- glu_panel(find_astro(seu,"GSE277487"),"TBI","Sham","GSE277487","CCI","重度","皮层","snRNA")
  rm(seu); gc(verbose=FALSE) }, error=function(e) cat("  GSE277487 失败:",conditionMessage(e),"\n"))

# ── GSE209552 人类重度 ──
cat("\n== GSE209552 (人类重度TBI/皮层/snRNA) ==\n"); flush.console()
tryCatch({ hum <- readRDS(HUM)
  hum$Group <- ifelse(grepl("HuBrainCTL",hum$sample_id),"Sham",
               ifelse(grepl("MJ_TBI",hum$sample_id),"TBI",NA))
  he <- GetAssayData(hum,layer="data",assay="RNA")
  mkh <- lapply(MK,toupper); cls <- sort(unique(hum$seurat_clusters))
  sc <- sapply(names(mkh),function(ct){ g <- intersect(mkh[[ct]],rownames(he))
    sapply(cls,function(cl) mean(as.matrix(he[g,colnames(hum)[hum$seurat_clusters==cl],drop=FALSE]))) })
  acl <- cls[colnames(sc)[apply(sc,1,which.max)]=="Astro"]
  hastro <- hum[,hum$seurat_clusters %in% acl]
  RES[["GSE209552"]] <- glu_panel(hastro,"TBI","Sham","GSE209552","TBI","人类重度","皮层","snRNA")
  rm(hum,hastro); gc(verbose=FALSE) }, error=function(e) cat("  GSE209552 失败:",conditionMessage(e),"\n"))

# ── GSE160763 FPI皮层（探索性混池）──
cat("\n== GSE160763 (FPI/皮层/scRNA·探索性) ==\n"); flush.console()
tryCatch({ d <- file.path(RAW,"GSE160763")
  gg <- as.list(c(GSM4878686="Sham",GSM4878687="Sham",GSM4878690="TBI",GSM4878691="TBI"))
  seu <- build_merged(d, gg)
  RES[["GSE160763"]] <- glu_panel(find_astro(seu,"GSE160763"),"TBI","Sham","GSE160763","FPI","探索·混池","皮层","scRNA")
  rm(seu); gc(verbose=FALSE) }, error=function(e) cat("  GSE160763 失败:",conditionMessage(e),"\n"))

# ── GSE302863 CCI皮层（探索性无重复）──
cat("\n== GSE302863 (CCI/皮层/snRNA·探索性) ==\n"); flush.console()
tryCatch({ d <- file.path(RAW,"GSE302863")
  gg <- as.list(c(GSM9112625="Sham",GSM9112626="TBI",GSM9112627="TBI",
                  GSM9112628="TBI",GSM9112629="TBI"))
  seu <- build_merged(d, gg)
  RES[["GSE302863"]] <- glu_panel(find_astro(seu,"GSE302863"),"TBI","Sham","GSE302863","CCI","探索·无重复","皮层","snRNA")
  rm(seu); gc(verbose=FALSE) }, error=function(e) cat("  GSE302863 失败:",conditionMessage(e),"\n"))

cat("\n\n══════════ 剩余 7 数据集结果 ══════════\n")
tab <- bind_rows(RES)
print(as.data.frame(tab), row.names=FALSE)
write.csv(tab, file.path(FIGOUT,"glutamate_panel_v2_part2.csv"), row.names=FALSE)
cat("\n========= 104b 完成 =========\n")
