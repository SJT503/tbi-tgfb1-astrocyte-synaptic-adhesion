# Script 104c: GSE271769 重度CCI皮层 谷氨酸面板（104b 因缺 hdf5r 跳过，补跑）
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(dplyr) })
options(future.globals.maxSize=4e9)
BASE <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project"
d <- file.path(BASE,"raw_data/GSE271769")
FIGOUT <- file.path(BASE,"wet_lab_data/results/figures/Fig_Homeostatic_Deepdive")

MK <- list(Neuron=c("Snap25","Syt1","Rbfox3"), Astro=c("Aqp4","Gfap","Slc1a2"),
  Oligo=c("Plp1","Mbp","Mobp"), Micro=c("Csf1r","P2ry12","Cx3cr1"),
  OPC=c("Pdgfra","Cspg4"), Endo=c("Cldn5","Flt1"))

find_astro <- function(seu) {
  if (ncol(seu) > 40000) { set.seed(1); seu <- seu[, sample(colnames(seu),40000)]
    cat(sprintf("    下采样到 40000\n")) }
  cat(sprintf("    聚类中 (n=%d)...\n", ncol(seu))); flush.console()
  seu <- NormalizeData(seu, verbose=FALSE)
  seu <- FindVariableFeatures(seu, nfeatures=2000, verbose=FALSE)
  seu <- ScaleData(seu, verbose=FALSE); seu <- RunPCA(seu, npcs=30, verbose=FALSE)
  seu <- FindNeighbors(seu, dims=1:30, verbose=FALSE)
  seu <- FindClusters(seu, resolution=0.5, verbose=FALSE)
  ex <- GetAssayData(seu, layer="data", assay="RNA")
  cls <- sort(unique(seu$seurat_clusters))
  sc <- sapply(names(MK), function(ct) {
    g <- intersect(MK[[ct]], rownames(ex))
    sapply(cls, function(cl) mean(as.matrix(ex[g, colnames(seu)[seu$seurat_clusters==cl], drop=FALSE]))) })
  acl <- cls[colnames(sc)[apply(sc,1,which.max)]=="Astro"]
  seu[, seu$seurat_clusters %in% acl]
}

cat("\n== GSE271769 (重度CCI/皮层/scRNA h5) ==\n"); flush.console()
h5 <- list.files(d, pattern="\\.h5$", full.names=TRUE)
objs <- lapply(h5, function(f) {
  grp <- if (grepl("NoTBI",basename(f))) "Sham" else "TBI"
  o <- CreateSeuratObject(Read10X_h5(f), min.cells=3, min.features=200)
  o$Group <- grp; o })
seu <- merge(objs[[1]],objs[-1]); seu <- JoinLayers(seu)
seu[["pm"]] <- PercentageFeatureSet(seu,pattern="^mt-")
seu <- subset(seu, nFeature_RNA>=200 & pm<25)
cat(sprintf("  QC后 %d 细胞 (Sham %d / TBI %d)\n",ncol(seu),sum(seu$Group=="Sham"),sum(seu$Group=="TBI")))
astro <- find_astro(seu)

grp <- astro$Group
counts <- GetAssayData(astro, layer="counts", assay="RNA")
med <- c(median(colSums(counts[,colnames(astro)[grp=="TBI"],drop=FALSE])),
         median(colSums(counts[,colnames(astro)[grp=="Sham"],drop=FALSE])))
set.seed(42)
ds <- SampleUMI(counts, max.umi=round(min(med)*0.9), upsample=FALSE, verbose=FALSE)
rownames(ds) <- rownames(counts); colnames(ds) <- colnames(counts)
dd <- CreateSeuratObject(counts=ds); dd$Group <- as.character(grp)
dd <- NormalizeData(dd, verbose=FALSE)
de <- GetAssayData(dd, layer="data", assay="RNA")
row <- data.frame(dataset="GSE271769", model="CCI", severity="重度", region="皮层",
                  platform="scRNA", n_astro=ncol(astro),
                  n_TBI=sum(grp=="TBI"), n_ctrl=sum(grp=="Sham"))
for (g in c("Glul","Slc1a2","Slc1a3")) {
  gg <- toupper(g)
  if (!g %in% rownames(de)) { row[[paste0(gg,"_pct")]] <- NA; row[[paste0(gg,"_p")]] <- NA; next }
  vt <- as.numeric(de[g, grp=="TBI"]); vc <- as.numeric(de[g, grp=="Sham"])
  row[[paste0(gg,"_pct")]] <- round((mean(vt)-mean(vc))/(mean(vc)+1e-9)*100,1)
  row[[paste0(gg,"_p")]] <- signif(tryCatch(wilcox.test(vt,vc)$p.value,error=function(e)NA),3)
}
cat(sprintf("\n  [GSE271769] 星形 %d (T%d/C%d) | GLUL %+.1f%% / SLC1A2 %+.1f%% / SLC1A3 %+.1f%%\n",
            ncol(astro), row$n_TBI, row$n_ctrl, row$GLUL_pct, row$SLC1A2_pct, row$SLC1A3_pct))
write.csv(row, file.path(FIGOUT,"glutamate_panel_v2_GSE271769.csv"), row.names=FALSE)
cat("\n========= 104c 完成 =========\n")
