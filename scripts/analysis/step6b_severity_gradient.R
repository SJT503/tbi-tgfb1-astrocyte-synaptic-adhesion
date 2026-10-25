#!/usr/bin/env Rscript
# Step 6b — neuronal survival / intrinsic-apoptosis / adhesion module scores
# across the TBI severity gradient: moderate (GSE290150), severe (GSE271769,
# GSE277487), human (GSE209552). Loaders reused from scripts/105.
# Real processed/raw data only.
suppressPackageStartupMessages({library(Seurat); library(Matrix); library(dplyr)})
options(future.globals.maxSize=4e9)
BASE <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project"
RAW  <- file.path(BASE,"raw_data")
HUM  <- file.path(BASE,"wet_lab_data/results/human_validation_v2/processed.rds")
OUT  <- file.path(BASE,"JNI/analysis/step6_out")
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

MODS <- list(
  IntrinsicApoptosis = c("Bax","Bak1","Bbc3","Pmaip1","Apaf1","Casp9","Casp3","Cycs"),
  Antiapoptotic      = c("Bcl2","Bcl2l1","Mcl1","Bcl2l2","Xiap"),
  Adhesion           = c("Ncam1","Nrxn1","Nrxn3","Cntn1","Negr1","Cadm1"))
MK <- list(Neuron=c("Snap25","Syt1","Rbfox3"), Astro=c("Aqp4","Gfap","Slc1a2"),
  Oligo=c("Plp1","Mbp","Mobp"), Micro=c("Csf1r","P2ry12","Cx3cr1"),
  OPC=c("Pdgfra","Cspg4"), Endo=c("Cldn5","Flt1"))

find_neuron <- function(seu, tag="") {
  if (ncol(seu) > 40000) { set.seed(1); seu <- seu[, sample(colnames(seu),40000)] }
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
    sapply(cls, function(cl) mean(as.matrix(ex[g, colnames(seu)[seu$seurat_clusters==cl],drop=FALSE]))) })
  ncl <- cls[colnames(sc)[apply(sc,1,which.max)]=="Neuron"]
  seu[, seu$seurat_clusters %in% ncl]
}

# depth-match TBI vs control, normalize, module score, test
score_neuron <- function(neu, dataset, severity) {
  upper <- toupper(rownames(neu)[1])==rownames(neu)[1]
  grp <- as.character(neu$Group); smp <- as.character(neu$sample)
  n1 <- sum(grp=="TBI",na.rm=TRUE); n0 <- sum(grp=="Sham",na.rm=TRUE)
  if (n1<20 || n0<20) { cat(sprintf("  [%s] neurons too few (%d/%d) skip\n",dataset,n1,n0)); return(NULL) }
  # purity gate
  ex0 <- GetAssayData(neu, layer="data", assay="RNA")
  sn <- if(upper)"SNAP25" else "Snap25"; aq <- if(upper)"AQP4" else "Aqp4"
  if (all(c(sn,aq) %in% rownames(ex0))) {
    r <- mean(ex0[sn,])/(mean(ex0[aq,])+1e-6)
    cat(sprintf("  [%s] PURITY Snap25/Aqp4=%.1f %s\n",dataset,r,ifelse(r>=3,"PASS","**CHECK**"))) }
  counts <- GetAssayData(neu, layer="counts", assay="RNA")
  med <- c(median(colSums(counts[,grp=="TBI",drop=FALSE])),
           median(colSums(counts[,grp=="Sham",drop=FALSE])))
  set.seed(42)
  ds <- SampleUMI(counts, max.umi=round(min(med)*0.9), upsample=FALSE, verbose=FALSE)
  rownames(ds) <- rownames(counts); colnames(ds) <- colnames(counts)
  d <- CreateSeuratObject(counts=ds); d$Group <- grp; d$sample <- smp
  d <- NormalizeData(d, verbose=FALSE)
  mods <- lapply(MODS, function(g){ gg <- if(upper) toupper(g) else g; gg[gg %in% rownames(d)] })
  set.seed(1); d <- AddModuleScore(d, features=mods, name="MOD_", assay="RNA")
  sc <- grep("^MOD_", colnames(d@meta.data), value=TRUE); names(sc) <- names(MODS)
  rows <- list()
  for (nm in names(MODS)) {
    v <- d@meta.data[[sc[nm]]]
    ct <- v[d$Group=="Sham"]; tb <- v[d$Group=="TBI"]
    # pseudobulk per sample
    pb <- tapply(v, d$sample, mean)
    pbg <- tapply(as.character(d$Group), d$sample, function(x) x[1])
    p0 <- pb[pbg=="Sham"]; p1 <- pb[pbg=="TBI"]
    pb_p <- if(length(p0)>=2 && length(p1)>=2) signif(t.test(p0,p1)$p.value,3) else NA
    rows[[nm]] <- data.frame(dataset=dataset, severity=severity, module=nm,
      n_ctl=length(ct), n_tbi=length(tb), n_samp_ctl=sum(pbg=="Sham"), n_samp_tbi=sum(pbg=="TBI"),
      mean_ctl=round(mean(ct),4), mean_tbi=round(mean(tb),4),
      delta=round(mean(tb)-mean(ct),4),
      percell_p=signif(wilcox.test(ct,tb)$p.value,3),
      pseudobulk_p=pb_p)
  }
  cat(sprintf("  [%s] neurons %d (T%d/C%d), samples C%d/T%d\n",
      dataset, ncol(neu), n1, n0, sum(tapply(as.character(d$Group),d$sample,function(x)x[1])=="Sham"),
      sum(tapply(as.character(d$Group),d$sample,function(x)x[1])=="TBI")))
  bind_rows(rows)
}

RES <- list()

cat("\n===== GSE290150 (moderate CCI cortex) =====\n"); flush.console()
tryCatch({ d <- file.path(RAW,"GSE290150")
  hdr <- scan(gzfile(file.path(d,"GSE290150_matrix.mtx.gz")), what=integer(), skip=1, nmax=3, quiet=TRUE)
  dt <- data.table::fread(file.path(d,"GSE290150_matrix.mtx.gz"), skip=2, header=FALSE)
  bc <- read.table(gzfile(file.path(d,"GSE290150_features.tsv.gz")),sep="\t")[,1]
  gn <- read.table(gzfile(file.path(d,"GSE290150_barcodes.tsv.gz")),sep="\t")[,1]
  m <- sparseMatrix(i=dt$V1,j=dt$V2,x=dt$V3,dims=hdr[1:2]); rm(dt); gc(verbose=FALSE)
  grpf <- sub("^[^.]*\\.","",bc); keep <- grepl("Normothermia",grpf)
  m <- t(m[keep,]); rownames(m) <- make.unique(as.character(gn)); colnames(m) <- bc[keep]
  gk <- grpf[keep]
  seu <- CreateSeuratObject(counts=m, min.cells=3, min.features=200)
  mm <- match(colnames(seu), bc[keep])
  seu$Group  <- ifelse(grepl("^Sham",gk[mm]),"Sham","TBI")
  seu$sample <- gk[mm]
  seu[["pm"]] <- PercentageFeatureSet(seu,pattern="^mt-")
  seu <- subset(seu, nFeature_RNA>=200 & pm<25)
  RES[["GSE290150"]] <- score_neuron(find_neuron(seu,"GSE290150"),"GSE290150","moderate")
  rm(m,seu); gc(verbose=FALSE) }, error=function(e) cat("GSE290150 FAIL:",conditionMessage(e),"\n"))

cat("\n===== GSE271769 (severe CCI cortex, h5) =====\n"); flush.console()
tryCatch({ d <- file.path(RAW,"GSE271769")
  h5 <- list.files(d, pattern="\\.h5$", full.names=TRUE)
  objs <- lapply(h5, function(f) {
    grp <- if (grepl("NoTBI",basename(f))) "Sham" else "TBI"
    o <- CreateSeuratObject(Read10X_h5(f), min.cells=3, min.features=200)
    o$Group <- grp; o$sample <- basename(f); o })
  seu <- merge(objs[[1]],objs[-1]); seu <- JoinLayers(seu)
  seu[["pm"]] <- PercentageFeatureSet(seu,pattern="^mt-")
  seu <- subset(seu, nFeature_RNA>=200 & pm<25)
  RES[["GSE271769"]] <- score_neuron(find_neuron(seu,"GSE271769"),"GSE271769","severe")
  rm(objs,seu); gc(verbose=FALSE) }, error=function(e) cat("GSE271769 FAIL:",conditionMessage(e),"\n"))

cat("\n===== GSE277487 (severe CCI cortex, snRNA) =====\n"); flush.console()
tryCatch({ d <- file.path(RAW,"GSE277487")
  gsms <- unique(sub("_.*$","",list.files(d,pattern="matrix\\.mtx")))
  objs <- lapply(gsms, function(x) {
    f <- list.files(d,pattern=paste0("^",x,"[_.]"),full.names=TRUE)
    o <- CreateSeuratObject(ReadMtx(mtx=grep("matrix",f,value=TRUE),cells=grep("barcodes",f,value=TRUE),
         features=grep("features|genes",f,value=TRUE),feature.column=2), min.cells=3,min.features=200)
    o$Group <- if (grepl("Ctrl",list.files(d,pattern=paste0("^",x,"_"))[1])) "Sham" else "TBI"
    o$sample <- x; o })
  seu <- merge(objs[[1]],objs[-1]); seu <- JoinLayers(seu)
  seu[["pm"]] <- PercentageFeatureSet(seu,pattern="^mt-")
  seu <- subset(seu, nFeature_RNA>=200 & pm<25)
  RES[["GSE277487"]] <- score_neuron(find_neuron(seu,"GSE277487"),"GSE277487","severe")
  rm(objs,seu); gc(verbose=FALSE) }, error=function(e) cat("GSE277487 FAIL:",conditionMessage(e),"\n"))

cat("\n===== GSE209552 (human severe TBI cortex) =====\n"); flush.console()
tryCatch({ hum <- readRDS(HUM)
  hum$Group <- ifelse(grepl("HuBrainCTL",hum$sample_id),"Sham",
               ifelse(grepl("MJ_TBI",hum$sample_id),"TBI",NA))
  hum$sample <- hum$sample_id
  hum <- hum[, !is.na(hum$Group)]
  he <- GetAssayData(hum,layer="data",assay="RNA")
  mkh <- lapply(MK,toupper); cls <- sort(unique(hum$seurat_clusters))
  sc <- sapply(names(mkh),function(ct){ g <- intersect(mkh[[ct]],rownames(he))
    sapply(cls,function(cl) mean(as.matrix(he[g,colnames(hum)[hum$seurat_clusters==cl],drop=FALSE]))) })
  ncl <- cls[colnames(sc)[apply(sc,1,which.max)]=="Neuron"]
  hneu <- hum[,hum$seurat_clusters %in% ncl]
  RES[["GSE209552"]] <- score_neuron(hneu,"GSE209552","human")
  rm(hum,hneu); gc(verbose=FALSE) }, error=function(e) cat("GSE209552 FAIL:",conditionMessage(e),"\n"))

tab <- bind_rows(RES)
write.csv(tab, file.path(OUT,"step6b_severity_gradient.csv"), row.names=FALSE)
cat("\n\n========== SEVERITY GRADIENT — neuronal module scores ==========\n")
print(as.data.frame(tab), row.names=FALSE)
cat("\nDONE step6b\n")
