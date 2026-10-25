# Script 105: PANoptosis 基因·神经元面板 — Step 1 免费计算门控
# ═══════════════════════════════════════════════════════════════════
# 目的：查 TBI 神经元是否上调 PANoptosis 三臂基因，决定要不要补焦亡湿实验
# 数据：自测35k图谱(描述性,批次混杂) + GSE290150/271769/277487(CCI皮层) + GSE209552(人源)
# ⚠️ 外部数据集无 TGFβRI-KO 臂；自测有3组但批次混杂→自测只作描述性

suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(dplyr) })
options(future.globals.maxSize=4e9)
BASE <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project"
RAW  <- file.path(BASE,"raw_data")
ATLAS <- file.path(BASE,"wet_lab_data/results/seurat/novogene_atlas_annotated.rds")
HUM  <- file.path(BASE,"wet_lab_data/results/human_validation_v2/processed.rds")
OUT  <- file.path(BASE,"wet_lab_data/results/figures/Fig_PANoptosis")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

# PANoptosis 三臂基因面板（小鼠符号；人源转大写）
PANO <- list(
  "坏死性凋亡" = c("Zbp1","Ripk1","Ripk3","Mlkl"),
  "凋亡"       = c("Casp8","Casp3","Casp9","Casp7","Casp6","Bax","Bcl2","Fadd"),
  "焦亡"       = c("Casp1","Nlrp3","Aim2","Pycard","Gsdmd","Gsdme","Il1b","Il18")
)
ALLG <- unlist(PANO, use.names=FALSE)
MK <- list(Neuron=c("Snap25","Syt1","Rbfox3"), Astro=c("Aqp4","Gfap","Slc1a2"),
  Oligo=c("Plp1","Mbp","Mobp"), Micro=c("Csf1r","P2ry12","Cx3cr1"),
  OPC=c("Pdgfra","Cspg4"), Endo=c("Cldn5","Flt1"))

find_neuron <- function(seu, tag="") {
  if (ncol(seu) > 40000) { set.seed(1); seu <- seu[, sample(colnames(seu),40000)]
    cat(sprintf("    [%s] 下采样到 40000\n",tag)) }
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
  ncl <- cls[colnames(sc)[apply(sc,1,which.max)]=="Neuron"]
  seu[, seu$seurat_clusters %in% ncl]
}

# 神经元 PANoptosis 面板：TBI vs 对照（深度匹配后）
pano_panel <- function(neu, g_tbi, g_ctrl, dataset) {
  upper <- toupper(rownames(neu)[1])==rownames(neu)[1]
  genes <- if (upper) toupper(ALLG) else ALLG
  grp <- neu$Group
  n1 <- sum(grp==g_tbi,na.rm=TRUE); n0 <- sum(grp==g_ctrl,na.rm=TRUE)
  if (n1<20 || n0<20) { cat(sprintf("  [%s] 神经元太少(%d/%d) 跳过\n",dataset,n1,n0)); return(NULL) }
  counts <- GetAssayData(neu, layer="counts", assay="RNA")
  med <- c(median(colSums(counts[,colnames(neu)[grp==g_tbi],drop=FALSE])),
           median(colSums(counts[,colnames(neu)[grp==g_ctrl],drop=FALSE])))
  set.seed(42)
  ds <- SampleUMI(counts, max.umi=round(min(med)*0.9), upsample=FALSE, verbose=FALSE)
  rownames(ds) <- rownames(counts); colnames(ds) <- colnames(counts)
  d <- CreateSeuratObject(counts=ds); d$Group <- as.character(grp)
  d <- NormalizeData(d, verbose=FALSE)
  de <- GetAssayData(d, layer="data", assay="RNA")
  res <- lapply(seq_along(ALLG), function(i) {
    g <- genes[i]
    if (!g %in% rownames(de)) return(data.frame(gene=ALLG[i],pct=NA,p=NA,pct_det_TBI=NA))
    vt <- as.numeric(de[g, grp==g_tbi]); vc <- as.numeric(de[g, grp==g_ctrl])
    data.frame(gene=ALLG[i],
               pct=round((mean(vt)-mean(vc))/(mean(vc)+1e-9)*100,1),
               p=signif(tryCatch(wilcox.test(vt,vc)$p.value,error=function(e)NA),3),
               pct_det_TBI=round(mean(vt>0)*100,1))
  })
  tab <- bind_rows(res)
  tab$category <- rep(names(PANO), sapply(PANO,length))
  tab$dataset <- dataset; tab$n_TBI <- n1; tab$n_ctrl <- n0
  cat(sprintf("  [%s] 神经元 %d (T%d/C%d)\n", dataset, ncol(neu), n1, n0))
  for (cn in names(PANO)) {
    sub <- tab[tab$category==cn,]
    up <- sum(sub$pct>0,na.rm=TRUE); sig <- sum(sub$p<0.05 & sub$pct>0,na.rm=TRUE)
    cat(sprintf("    %s: %d/%d 基因↑，其中显著↑ %d\n", cn, up, sum(!is.na(sub$pct)), sig))
  }
  tab
}

RES <- list()

# ═══ 自测 35k 图谱（描述性，3 组）═══
cat("\n========== 自测 35k 图谱（描述性·批次混杂）==========\n"); flush.console()
tryCatch({ s <- readRDS(ATLAS)
  neu <- s[, s$novo_broad=="Neuron"]
  ex <- GetAssayData(neu, layer="data", assay="RNA")
  cat(sprintf("自测神经元 %d (Sham %d / CCI_Scr %d / KO %d)\n", ncol(neu),
      sum(neu$stim=="Sham"), sum(neu$stim=="CCI_Scr"), sum(neu$stim=="CCI_TGFbetaR1_KO")))
  selfdesc <- lapply(ALLG, function(g) {
    if (!g %in% rownames(ex)) return(data.frame(gene=g,Sham=NA,CCI_Scr=NA,CCI_KO=NA,det_CCI=NA))
    data.frame(gene=g,
      Sham=round(mean(ex[g,neu$stim=="Sham"]),3),
      CCI_Scr=round(mean(ex[g,neu$stim=="CCI_Scr"]),3),
      CCI_KO=round(mean(ex[g,neu$stim=="CCI_TGFbetaR1_KO"]),3),
      det_CCI=round(mean(ex[g,neu$stim=="CCI_Scr"]>0)*100,1))
  })
  sd_tab <- bind_rows(selfdesc)
  sd_tab$category <- rep(names(PANO), sapply(PANO,length))
  cat("\n自测神经元 PANoptosis 基因（per-group 均值，⚠️批次混杂→仅描述性，不作定量）：\n")
  print(sd_tab[,c("category","gene","Sham","CCI_Scr","CCI_KO","det_CCI")], row.names=FALSE)
  write.csv(sd_tab, file.path(OUT,"panoptosis_self_atlas_descriptive.csv"), row.names=FALSE)
  rm(s,neu); gc(verbose=FALSE) }, error=function(e) cat("自测图谱失败:",conditionMessage(e),"\n"))

# ═══ GSE290150 中度CCI ═══
cat("\n========== GSE290150（中度CCI皮层）==========\n"); flush.console()
tryCatch({ d <- file.path(RAW,"GSE290150")
  hdr <- scan(gzfile(file.path(d,"GSE290150_matrix.mtx.gz")), what=integer(), skip=1, nmax=3, quiet=TRUE)
  dt <- data.table::fread(file.path(d,"GSE290150_matrix.mtx.gz"), skip=2, header=FALSE)
  bc <- read.table(gzfile(file.path(d,"GSE290150_features.tsv.gz")),sep="\t")[,1]
  gn <- read.table(gzfile(file.path(d,"GSE290150_barcodes.tsv.gz")),sep="\t")[,1]
  m <- sparseMatrix(i=dt$V1,j=dt$V2,x=dt$V3,dims=hdr[1:2]); rm(dt); gc(verbose=FALSE)
  grpf <- sub("^[^.]*\\.","",bc); keep <- grepl("Normothermia",grpf)
  m <- t(m[keep,]); rownames(m) <- make.unique(as.character(gn)); colnames(m) <- bc[keep]
  seu <- CreateSeuratObject(counts=m, min.cells=3, min.features=200)
  seu$Group <- ifelse(grepl("^Sham",grpf[keep][match(colnames(seu),bc[keep])]),"Sham","TBI")
  seu[["pm"]] <- PercentageFeatureSet(seu,pattern="^mt-")
  seu <- subset(seu, nFeature_RNA>=200 & pm<25)
  RES[["GSE290150"]] <- pano_panel(find_neuron(seu,"GSE290150"),"TBI","Sham","GSE290150")
  rm(m,seu); gc(verbose=FALSE) }, error=function(e) cat("GSE290150 失败:",conditionMessage(e),"\n"))

# ═══ GSE271769 重度CCI (h5) ═══
cat("\n========== GSE271769（重度CCI皮层）==========\n"); flush.console()
tryCatch({ d <- file.path(RAW,"GSE271769")
  h5 <- list.files(d, pattern="\\.h5$", full.names=TRUE)
  objs <- lapply(h5, function(f) {
    grp <- if (grepl("NoTBI",basename(f))) "Sham" else "TBI"
    o <- CreateSeuratObject(Read10X_h5(f), min.cells=3, min.features=200); o$Group <- grp; o })
  seu <- merge(objs[[1]],objs[-1]); seu <- JoinLayers(seu)
  seu[["pm"]] <- PercentageFeatureSet(seu,pattern="^mt-")
  seu <- subset(seu, nFeature_RNA>=200 & pm<25)
  RES[["GSE271769"]] <- pano_panel(find_neuron(seu,"GSE271769"),"TBI","Sham","GSE271769")
  rm(objs,seu); gc(verbose=FALSE) }, error=function(e) cat("GSE271769 失败:",conditionMessage(e),"\n"))

# ═══ GSE277487 重度CCI snRNA ═══
cat("\n========== GSE277487（重度CCI皮层·snRNA）==========\n"); flush.console()
tryCatch({ d <- file.path(RAW,"GSE277487")
  gsms <- unique(sub("_.*$","",list.files(d,pattern="matrix\\.mtx")))
  objs <- lapply(gsms, function(x) {
    f <- list.files(d,pattern=paste0("^",x,"[_.]"),full.names=TRUE)
    o <- CreateSeuratObject(ReadMtx(mtx=grep("matrix",f,value=TRUE),cells=grep("barcodes",f,value=TRUE),
         features=grep("features|genes",f,value=TRUE),feature.column=2), min.cells=3,min.features=200)
    o$Group <- if (grepl("Ctrl",list.files(d,pattern=paste0("^",x,"_"))[1])) "Sham" else "TBI"; o })
  seu <- merge(objs[[1]],objs[-1]); seu <- JoinLayers(seu)
  seu[["pm"]] <- PercentageFeatureSet(seu,pattern="^mt-")
  seu <- subset(seu, nFeature_RNA>=200 & pm<25)
  RES[["GSE277487"]] <- pano_panel(find_neuron(seu,"GSE277487"),"TBI","Sham","GSE277487")
  rm(objs,seu); gc(verbose=FALSE) }, error=function(e) cat("GSE277487 失败:",conditionMessage(e),"\n"))

# ═══ GSE209552 人类重度 ═══
cat("\n========== GSE209552（人类重度TBI皮层）==========\n"); flush.console()
tryCatch({ hum <- readRDS(HUM)
  hum$Group <- ifelse(grepl("HuBrainCTL",hum$sample_id),"Sham",
               ifelse(grepl("MJ_TBI",hum$sample_id),"TBI",NA))
  he <- GetAssayData(hum,layer="data",assay="RNA")
  mkh <- lapply(MK,toupper); cls <- sort(unique(hum$seurat_clusters))
  sc <- sapply(names(mkh),function(ct){ g <- intersect(mkh[[ct]],rownames(he))
    sapply(cls,function(cl) mean(as.matrix(he[g,colnames(hum)[hum$seurat_clusters==cl],drop=FALSE]))) })
  ncl <- cls[colnames(sc)[apply(sc,1,which.max)]=="Neuron"]
  hneu <- hum[,hum$seurat_clusters %in% ncl]
  RES[["GSE209552"]] <- pano_panel(hneu,"TBI","Sham","GSE209552")
  rm(hum,hneu); gc(verbose=FALSE) }, error=function(e) cat("GSE209552 失败:",conditionMessage(e),"\n"))

# ═══ 汇总 ═══
cat("\n\n══════════ 外部数据集 神经元 PANoptosis 面板（TBI vs 对照）══════════\n")
tab <- bind_rows(RES)
if (nrow(tab)>0) {
  print(as.data.frame(tab[,c("dataset","category","gene","pct","p","pct_det_TBI")]), row.names=FALSE)
  write.csv(tab, file.path(OUT,"panoptosis_external_neuron_panel.csv"), row.names=FALSE)
}
cat("\n判读：三臂(坏死性凋亡/凋亡/焦亡)若都在 TBI 神经元↑ → 支持 PANoptosis 方向，焦亡湿实验值得补\n")
cat("\n========= 105 完成 =========\n")
