# Script 106: Kir4.1(Kcnj10)/AQP4 跨数据集门控 + Smad 活性相关
# ═══════════════════════════════════════════════════════════════════
# 目的：免费计算门控——TGFβRI→Kir4.1 机制候选
#   ① 星形 Kcnj10/Aqp4 在 TBI 是否下调（跨损伤梯度）
#   ② Smad 靶基因签名活性 vs Kcnj10 是否负相关（支持"TGFβ 信号压制 Kir4.1"）
#   ③ 自测图谱描述性：Kcnj10 在 Sham/CCI_Scr/CCI_KO 三组（KO 后是否回升）
# ⚠️ 自测 snRNA 批次混杂→自测只作描述性；外部无 KO 臂
# 数据：自测图谱 + 轻/中/重 CCI/FPI/TBI + 人源

suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(dplyr) })
options(future.globals.maxSize=4e9)
BASE <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project"
RAW  <- file.path(BASE,"raw_data")
ATLAS <- file.path(BASE,"wet_lab_data/results/seurat/novogene_atlas_annotated.rds")
HUM  <- file.path(BASE,"wet_lab_data/results/human_validation_v2/processed.rds")
OUT  <- file.path(BASE,"wet_lab_data/results/figures/Fig_Kir41_Gate")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

# Smad/TGFβ 经典转录靶基因（同 Script 102）
SMAD_MM <- c("Serpine1","Smad7","Skil","Pmepa1","Junb","Id1","Id2","Id3",
             "Tgfbi","Ltbp2","Cdkn1a","Smad6","Gadd45b","Bcl3")
MK <- list(Neuron=c("Snap25","Syt1","Rbfox3"), Astro=c("Aqp4","Gfap","Slc1a2"),
  Oligo=c("Plp1","Mbp","Mobp"), Micro=c("Csf1r","P2ry12","Cx3cr1"),
  OPC=c("Pdgfra","Cspg4"), Endo=c("Cldn5","Flt1"))

find_astro <- function(seu, tag="") {
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
  acl <- cls[colnames(sc)[apply(sc,1,which.max)]=="Astro"]
  seu[, seu$seurat_clusters %in% acl]
}

# Kir4.1/AQP4 面板 + Smad 相关
kir_panel <- function(astro, g_tbi, g_ctrl, dataset, severity) {
  upper <- toupper(rownames(astro)[1])==rownames(astro)[1]
  KIR <- if (upper) "KCNJ10" else "Kcnj10"
  AQP <- if (upper) "AQP4" else "Aqp4"
  smad <- if (upper) toupper(SMAD_MM) else SMAD_MM
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
  row <- data.frame(dataset=dataset, severity=severity, n_astro=ncol(astro), n_TBI=n1, n_ctrl=n0)
  for (g in c(KIR,AQP)) {
    nm <- if (g %in% c("Kcnj10","KCNJ10")) "Kir41" else "Aqp4"
    if (!g %in% rownames(de)) { row[[paste0(nm,"_pct")]] <- NA; row[[paste0(nm,"_p")]] <- NA
      row[[paste0(nm,"_det")]] <- NA; next }
    vt <- as.numeric(de[g, grp==g_tbi]); vc <- as.numeric(de[g, grp==g_ctrl])
    row[[paste0(nm,"_pct")]] <- round((mean(vt)-mean(vc))/(mean(vc)+1e-9)*100,1)
    row[[paste0(nm,"_p")]]   <- signif(tryCatch(wilcox.test(vt,vc)$p.value,error=function(e)NA),3)
    row[[paste0(nm,"_det")]] <- round(mean(vt>0)*100,1)
  }
  # Smad 活性 vs Kcnj10 相关（全部星形）
  sg <- intersect(smad, rownames(de))
  if (length(sg)>=5 && KIR %in% rownames(de)) {
    d <- AddModuleScore(d, features=list(sg), name="Smad", assay="RNA")
    sc <- d$Smad1; kv <- as.numeric(de[KIR,])
    ct <- suppressWarnings(cor.test(sc, kv, method="spearman"))
    row$Smad_Kir_r <- round(ct$estimate,3)
    row$Smad_Kir_p <- signif(ct$p.value,3)
    row$smad_hits  <- length(sg)
  } else { row$Smad_Kir_r <- NA; row$Smad_Kir_p <- NA; row$smad_hits <- length(sg) }
  cat(sprintf("  [%s] 星形%d (T%d/C%d) | Kir4.1 %+.1f%%(p=%.1e,det%.1f%%) | AQP4 %+.1f%% | Smad-Kir r=%+.3f\n",
      dataset, ncol(astro), n1, n0, row$Kir41_pct, row$Kir41_p, row$Kir41_det,
      row$Aqp4_pct, ifelse(is.na(row$Smad_Kir_r),0,row$Smad_Kir_r)))
  row
}

build_merged <- function(d, gsm_grp, mt="^mt-") {
  objs <- lapply(names(gsm_grp), function(gsm) {
    f <- list.files(d, pattern=paste0("^",gsm,"[_.]"), full.names=TRUE)
    mtx <- grep("matrix\\.mtx",f,value=TRUE); bc <- grep("barcodes\\.tsv",f,value=TRUE)
    ft <- grep("(features|genes)\\.tsv",f,value=TRUE)
    nc <- length(strsplit(readLines(gzfile(ft),n=1),"\t")[[1]])
    o <- CreateSeuratObject(ReadMtx(mtx=mtx,cells=bc,features=ft,feature.column=min(nc,2)),
         min.cells=3,min.features=200)
    o$Group <- gsm_grp[[gsm]]; o })
  seu <- if (length(objs)==1) objs[[1]] else merge(objs[[1]],objs[-1])
  if (length(objs)>1) seu <- JoinLayers(seu)
  seu[["pm"]] <- PercentageFeatureSet(seu, pattern=mt)
  subset(seu, nFeature_RNA>=200 & pm<25)
}

RES <- list()

# ═══ 自测 35k 图谱（描述性，含 KO 臂）═══
cat("\n========== 自测 35k 图谱（描述性·含 KO 臂）==========\n"); flush.console()
tryCatch({ s <- readRDS(ATLAS)
  astro <- s[, s$novo_broad=="Astrocytes"]
  ex <- GetAssayData(astro, layer="data", assay="RNA")
  cat(sprintf("自测星形 %d (Sham %d/CCI_Scr %d/KO %d)\n", ncol(astro),
      sum(astro$stim=="Sham"),sum(astro$stim=="CCI_Scr"),sum(astro$stim=="CCI_TGFbetaR1_KO")))
  for (g in c("Kcnj10","Aqp4")) {
    if (!g %in% rownames(ex)) { cat(sprintf("  %s 不在\n",g)); next }
    cat(sprintf("  %s | Sham %.3f | CCI_Scr %.3f | CCI_KO %.3f | det_CCI %.1f%%\n", g,
        mean(ex[g,astro$stim=="Sham"]), mean(ex[g,astro$stim=="CCI_Scr"]),
        mean(ex[g,astro$stim=="CCI_TGFbetaR1_KO"]),
        mean(ex[g,astro$stim=="CCI_Scr"]>0)*100))
  }
  # Smad-Kir 相关（CCI_Scr 星形）
  sg <- intersect(SMAD_MM, rownames(ex))
  cci <- astro[, astro$stim=="CCI_Scr"]
  cci <- AddModuleScore(cci, features=list(sg), name="Smad", assay="RNA")
  if ("Kcnj10" %in% rownames(ex)) {
    kv <- as.numeric(GetAssayData(cci,layer="data",assay="RNA")["Kcnj10",])
    ct <- suppressWarnings(cor.test(cci$Smad1, kv, method="spearman"))
    cat(sprintf("  自测 CCI_Scr 星形：Smad活性 vs Kcnj10  r=%+.3f (p=%.1e)\n", ct$estimate, ct$p.value))
  }
  rm(s,astro); gc(verbose=FALSE) }, error=function(e) cat("自测图谱失败:",conditionMessage(e),"\n"))

# ═══ 人源 GSE209552（先跑，快）═══
cat("\n========== GSE209552（人类重度TBI）==========\n"); flush.console()
tryCatch({ hum <- readRDS(HUM)
  hum$Group <- ifelse(grepl("HuBrainCTL",hum$sample_id),"Sham",
               ifelse(grepl("MJ_TBI",hum$sample_id),"TBI",NA))
  he <- GetAssayData(hum,layer="data",assay="RNA")
  mkh <- lapply(MK,toupper); cls <- sort(unique(hum$seurat_clusters))
  sc <- sapply(names(mkh),function(ct){ g <- intersect(mkh[[ct]],rownames(he))
    sapply(cls,function(cl) mean(as.matrix(he[g,colnames(hum)[hum$seurat_clusters==cl],drop=FALSE]))) })
  acl <- cls[colnames(sc)[apply(sc,1,which.max)]=="Astro"]
  hastro <- hum[,hum$seurat_clusters %in% acl]
  RES[["GSE209552"]] <- kir_panel(hastro,"TBI","Sham","GSE209552","人类重度")
  rm(hum,hastro); gc(verbose=FALSE) }, error=function(e) cat("GSE209552 失败:",conditionMessage(e),"\n"))

# ═══ GSE101901 轻度TBI海马 ═══
cat("\n========== GSE101901（轻度TBI/海马）==========\n"); flush.console()
tryCatch({ mat <- as.matrix(read.table(gzfile(file.path(RAW,"GSE101901/GSE101901_DropSeqTBI.digital_expression.txt.gz")),
    header=TRUE,row.names=1,sep="\t",check.names=FALSE))
  seu <- CreateSeuratObject(counts=as(mat,"dgCMatrix"),min.cells=3,min.features=200)
  seu$Group <- ifelse(grepl("^Sham",sub("_.*$","",colnames(seu))),"Sham","TBI")
  seu[["pm"]] <- PercentageFeatureSet(seu,pattern="^mt-")
  seu <- subset(seu, nFeature_RNA>=200 & pm<25)
  RES[["GSE101901"]] <- kir_panel(find_astro(seu,"GSE101901"),"TBI","Sham","GSE101901","轻度")
  rm(mat,seu); gc(verbose=FALSE) }, error=function(e) cat("GSE101901 失败:",conditionMessage(e),"\n"))

# ═══ GSE180862 轻度FPI 皮层+海马 ═══
for (rg in c("Cortex","Hippocampus")) {
  cat(sprintf("\n========== GSE180862_%s（轻度FPI）==========\n",rg)); flush.console()
  tryCatch({ pre <- file.path(RAW,"GSE180862",paste0("GSE180862_DropSeq.",rg,".7days.24hrs.TBI.Sham"))
  mat <- readMM(gzfile(paste0(pre,".digital_expression.mtx.gz")))
  ft <- read.table(gzfile(paste0(pre,".features.tsv.gz")),sep="\t")[,1]
  bc <- read.table(gzfile(paste0(pre,".barcodes.tsv.gz")),sep="\t")[,1]
  if (nrow(mat)!=length(ft)) mat <- t(mat)
  rownames(mat) <- ft; colnames(mat) <- bc
  meta <- read.table(gzfile(paste0(pre,".metaData.tsv.gz")),header=TRUE,row.names=1,sep="\t")
  cm <- intersect(colnames(mat),rownames(meta)); mat <- mat[,cm]; meta <- meta[cm,]
  seu <- CreateSeuratObject(counts=as(mat,"dgCMatrix"),meta.data=meta)
  astro <- seu[,seu$CellType=="ASC"]; astro$Group <- astro$Condition
  RES[[paste0("GSE180862_",substr(rg,1,3))]] <- kir_panel(astro,"TBI","Sham",
    paste0("GSE180862_",substr(rg,1,3)),"轻度")
  rm(mat,seu,astro,meta); gc(verbose=FALSE) }, error=function(e) cat("失败:",conditionMessage(e),"\n"))
}

# ═══ GSE290150 中度CCI皮层 ═══
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
  RES[["GSE290150"]] <- kir_panel(find_astro(seu,"GSE290150"),"TBI","Sham","GSE290150","中度CCI")
  rm(m,seu); gc(verbose=FALSE) }, error=function(e) cat("GSE290150 失败:",conditionMessage(e),"\n"))

# ═══ GSE271769 重度CCI皮层 (h5) ═══
cat("\n========== GSE271769（重度CCI皮层）==========\n"); flush.console()
tryCatch({ d <- file.path(RAW,"GSE271769")
  h5 <- list.files(d, pattern="\\.h5$", full.names=TRUE)
  objs <- lapply(h5, function(f) {
    grp <- if (grepl("NoTBI",basename(f))) "Sham" else "TBI"
    o <- CreateSeuratObject(Read10X_h5(f), min.cells=3, min.features=200); o$Group <- grp; o })
  seu <- merge(objs[[1]],objs[-1]); seu <- JoinLayers(seu)
  seu[["pm"]] <- PercentageFeatureSet(seu,pattern="^mt-")
  seu <- subset(seu, nFeature_RNA>=200 & pm<25)
  RES[["GSE271769"]] <- kir_panel(find_astro(seu,"GSE271769"),"TBI","Sham","GSE271769","重度CCI")
  rm(objs,seu); gc(verbose=FALSE) }, error=function(e) cat("GSE271769 失败:",conditionMessage(e),"\n"))

# ═══ GSE277487 重度CCI皮层 (snRNA) ═══
cat("\n========== GSE277487（重度CCI皮层·snRNA）==========\n"); flush.console()
tryCatch({ d <- file.path(RAW,"GSE277487")
  gsms <- unique(sub("_.*$","",list.files(d,pattern="matrix\\.mtx")))
  gg <- setNames(lapply(gsms, function(x) {
    if (grepl("Ctrl",list.files(d,pattern=paste0("^",x,"_"))[1])) "Sham" else "TBI" }), gsms)
  seu <- build_merged(d, gg)
  RES[["GSE277487"]] <- kir_panel(find_astro(seu,"GSE277487"),"TBI","Sham","GSE277487","重度CCI")
  rm(seu); gc(verbose=FALSE) }, error=function(e) cat("GSE277487 失败:",conditionMessage(e),"\n"))

# ═══ 汇总 ═══
cat("\n\n══════════ Kir4.1/AQP4 跨数据集门控汇总 ══════════\n")
tab <- bind_rows(RES)
if (nrow(tab)>0) { print(as.data.frame(tab), row.names=FALSE)
  write.csv(tab, file.path(OUT,"kir41_crossdataset_gate.csv"), row.names=FALSE) }
cat("\n判读：\n")
cat(" ① Kir4.1(Kcnj10) 若跨数据集一致↓ → 门控通过第一关\n")
cat(" ② Smad活性 vs Kcnj10 若负相关(r<0) → 支持'TGFβ信号压制Kir4.1'方向\n")
cat(" ③ 自测：CCI_KO 的 Kcnj10 若比 CCI_Scr 回升 → 提示性支持(批次混杂,非定论)\n")
cat("\n========= 106 完成 =========\n")
