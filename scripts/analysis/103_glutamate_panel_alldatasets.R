# Script 103: 全数据集星形谷氨酸处理面板 — GLT-1/GLAST/GS 跨损伤模型对比
# ═══════════════════════════════════════════════════════════════════
# 统一方法（marker定星形 / ASC标注、深度降采样拉平、TBI-vs-对照）
# 数据集：轻度TBI、轻度FPI(皮层+海马)、轻度CCI、爆震、中度CCI、人类
# 只看 3 个谷氨酸处理基因：GLUL(GS) / Slc1a2(GLT-1) / Slc1a3(GLAST)

suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(dplyr) })

BASE <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project"
RAW  <- file.path(BASE, "raw_data")
G226 <- file.path(BASE, "wet_lab_data/external_data/GSE226072")
HUM  <- file.path(BASE, "wet_lab_data/results/human_validation_v2/processed.rds")
FIGOUT <- file.path(BASE, "wet_lab_data/results/figures/Fig_Homeostatic_Deepdive")

MK <- list(Neuron=c("Snap25","Syt1","Rbfox3"), Astro=c("Aqp4","Gfap","Slc1a2"),
  Oligo=c("Plp1","Mbp","Mobp"), Micro=c("Csf1r","P2ry12","Cx3cr1"),
  OPC=c("Pdgfra","Cspg4"), Endo=c("Cldn5","Flt1"))

find_astro <- function(seu) {
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

# 谷氨酸面板：深度降采样后算 3 基因 TBI-vs-对照
glu_panel <- function(astro, g_tbi, g_ctrl, dataset, model, region, platform) {
  upper <- toupper(rownames(astro)[1])==rownames(astro)[1]
  GLU <- if (upper) c("GLUL","SLC1A2","SLC1A3") else c("Glul","Slc1a2","Slc1a3")
  grp <- astro$Group
  n1 <- sum(grp==g_tbi); n0 <- sum(grp==g_ctrl)
  if (n1<15 || n0<15) { cat(sprintf("  [%s] 星形太少 跳过\n",dataset)); return(NULL) }
  counts <- GetAssayData(astro, layer="counts", assay="RNA")
  med <- c(median(colSums(counts[,colnames(astro)[grp==g_tbi],drop=FALSE])),
           median(colSums(counts[,colnames(astro)[grp==g_ctrl],drop=FALSE])))
  set.seed(42)
  ds <- SampleUMI(counts, max.umi=round(min(med)*0.9), upsample=FALSE, verbose=FALSE)
  rownames(ds) <- rownames(counts); colnames(ds) <- colnames(counts)
  d <- CreateSeuratObject(counts=ds); d$Group <- as.character(grp)
  d <- NormalizeData(d, verbose=FALSE)
  de <- GetAssayData(d, layer="data", assay="RNA")
  row <- data.frame(dataset=dataset, model=model, region=region, platform=platform,
                    n_astro=ncol(astro), n_TBI=n1, n_ctrl=n0)
  for (g in GLU) {
    gg <- toupper(g)
    if (!g %in% rownames(de)) { row[[paste0(gg,"_pct")]] <- NA; row[[paste0(gg,"_p")]] <- NA; next }
    vt <- as.numeric(de[g, grp==g_tbi]); vc <- as.numeric(de[g, grp==g_ctrl])
    row[[paste0(gg,"_pct")]] <- round((mean(vt)-mean(vc))/(mean(vc)+1e-9)*100,1)
    row[[paste0(gg,"_p")]] <- signif(tryCatch(wilcox.test(vt,vc)$p.value,error=function(e)NA),3)
  }
  cat(sprintf("  [%s] 星形 %d | GLUL %+.1f%% / SLC1A2 %+.1f%% / SLC1A3 %+.1f%%\n",
              dataset, ncol(astro), row$GLUL_pct, row$SLC1A2_pct, row$SLC1A3_pct))
  row
}

RES <- list()

# ── GSE101901 轻度TBI海马 ──
cat("\n== GSE101901 (轻度TBI/海马/DropSeq) ==\n")
{ mat <- as.matrix(read.table(gzfile(file.path(RAW,"GSE101901/GSE101901_DropSeqTBI.digital_expression.txt.gz")),
    header=TRUE,row.names=1,sep="\t",check.names=FALSE))
  seu <- CreateSeuratObject(counts=as(mat,"dgCMatrix"),min.cells=3,min.features=200)
  seu$Group <- ifelse(grepl("^Sham",sub("_.*$","",colnames(seu))),"Sham","TBI")
  seu[["pm"]] <- PercentageFeatureSet(seu,pattern="^mt-")
  seu <- subset(seu, nFeature_RNA>=200 & pm<25)
  RES[["GSE101901"]] <- glu_panel(find_astro(seu),"TBI","Sham","GSE101901","轻度TBI","海马","DropSeq")
  rm(mat,seu); gc(verbose=FALSE) }

# ── GSE180862 轻度FPI 皮层+海马 ──
for (rg in c("Cortex","Hippocampus")) {
  cat(sprintf("\n== GSE180862_%s (轻度FPI/DropSeq) ==\n",rg))
  pre <- file.path(RAW,"GSE180862",paste0("GSE180862_DropSeq.",rg,".7days.24hrs.TBI.Sham"))
  mat <- readMM(gzfile(paste0(pre,".digital_expression.mtx.gz")))
  ft <- read.table(gzfile(paste0(pre,".features.tsv.gz")),sep="\t")[,1]
  bc <- read.table(gzfile(paste0(pre,".barcodes.tsv.gz")),sep="\t")[,1]
  if (nrow(mat)!=length(ft)) mat <- t(mat)
  rownames(mat) <- ft; colnames(mat) <- bc
  meta <- read.table(gzfile(paste0(pre,".metaData.tsv.gz")),header=TRUE,row.names=1,sep="\t")
  cm <- intersect(colnames(mat),rownames(meta)); mat <- mat[,cm]; meta <- meta[cm,]
  seu <- CreateSeuratObject(counts=as(mat,"dgCMatrix"),meta.data=meta)
  astro <- seu[,seu$CellType=="ASC"]; astro$Group <- astro$Condition
  RES[[paste0("GSE180862_",substr(rg,1,3))]] <- glu_panel(astro,"TBI","Sham",
    paste0("GSE180862_",substr(rg,1,3)),"轻度FPI",rg,"DropSeq")
  rm(mat,seu,astro,meta); gc(verbose=FALSE)
}

# ── GSE292715 轻度CCI海马 / GSE230253 爆震海马（1v1）──
for (dd in list(c("GSE292715","轻度CCI","GSM8864084_c-HIP","GSM8864083_t-HIP"),
                c("GSE230253","爆震TBI","GSM7210822_matrix_c-TBI|GSM7210822_barcodes_c-TBI|GSM7210822_genes_c-TBI",
                  "GSM7210821_matrix_b-TBI|GSM7210821_barcodes_b-TBI|GSM7210821_genes_b-TBI"))) {
  cat(sprintf("\n== %s (%s/海马/10x) ==\n",dd[1],dd[2]))
  rdmtx <- function(pfx) {
    if (grepl("\\|",pfx)) { p <- strsplit(pfx,"\\|")[[1]]
      ReadMtx(mtx=file.path(RAW,dd[1],paste0(p[1],".mtx.gz")),
              cells=file.path(RAW,dd[1],paste0(p[2],".tsv.gz")),
              features=file.path(RAW,dd[1],paste0(p[3],".tsv.gz")),feature.column=2)
    } else ReadMtx(mtx=file.path(RAW,dd[1],paste0(pfx,"_matrix.mtx.gz")),
              cells=file.path(RAW,dd[1],paste0(pfx,"_barcodes.tsv.gz")),
              features=file.path(RAW,dd[1],paste0(pfx,"_genes.tsv.gz")),feature.column=2) }
  oc <- CreateSeuratObject(rdmtx(dd[3]),min.cells=3,min.features=200); oc$Group <- "Sham"
  ot <- CreateSeuratObject(rdmtx(dd[4]),min.cells=3,min.features=200); ot$Group <- "TBI"
  seu <- merge(oc,ot); seu <- JoinLayers(seu)
  seu[["pm"]] <- PercentageFeatureSet(seu,pattern="^mt-")
  seu <- subset(seu, nFeature_RNA>=200 & pm<25)
  RES[[dd[1]]] <- glu_panel(find_astro(seu),"TBI","Sham",dd[1],dd[2],"海马","10x")
  rm(oc,ot,seu); gc(verbose=FALSE)
}

# ── GSE226072 中度CCI皮层 ──
cat("\n== GSE226072 (中度CCI/皮层/snRNA) ==\n")
{ S <- list(c("GSM7062619_D1_S1","Sham"),c("GSM7062620_D1_S2","CCI"),
            c("GSM7062623_D1_S5","Sham"),c("GSM7062624_D1_S6","CCI"),
            c("GSM7062627_D1_S9","Sham"),c("GSM7062628_D1_S10","CCI"),
            c("GSM7062631_D1_S13","Sham"),c("GSM7062632_D1_S14","CCI"))
  objs <- lapply(S,function(s){ m <- ReadMtx(mtx=file.path(G226,paste0(s[1],"matrix.mtx.gz")),
    cells=file.path(G226,paste0(s[1],"barcodes.tsv.gz")),
    features=file.path(G226,paste0(s[1],"features.tsv.gz")),feature.column=2)
    o <- CreateSeuratObject(m,min.cells=3,min.features=200); o$Group <- s[2]; o })
  seu <- merge(objs[[1]],objs[-1]); seu <- JoinLayers(seu)
  seu[["pm"]] <- PercentageFeatureSet(seu,pattern="^mt-")
  seu <- subset(seu, nFeature_RNA>=200 & pm<25)
  RES[["GSE226072"]] <- glu_panel(find_astro(seu),"CCI","Sham","GSE226072","中度CCI","损伤旁皮层","snRNA")
  rm(objs,seu); gc(verbose=FALSE) }

# ── 人类 GSE209552 ──
cat("\n== GSE209552 (人类重度TBI/皮层/snRNA) ==\n")
{ hum <- readRDS(HUM)
  hum$Group <- ifelse(grepl("HuBrainCTL",hum$sample_id),"Sham",
               ifelse(grepl("MJ_TBI",hum$sample_id),"TBI",NA))
  he <- GetAssayData(hum,layer="data",assay="RNA")
  mkh <- lapply(MK,toupper); cls <- sort(unique(hum$seurat_clusters))
  sc <- sapply(names(mkh),function(ct){ g <- intersect(mkh[[ct]],rownames(he))
    sapply(cls,function(cl) mean(as.matrix(he[g,colnames(hum)[hum$seurat_clusters==cl],drop=FALSE]))) })
  acl <- cls[colnames(sc)[apply(sc,1,which.max)]=="Astro"]
  hastro <- hum[,hum$seurat_clusters %in% acl]
  RES[["GSE209552"]] <- glu_panel(hastro,"TBI","Sham","GSE209552","人类重度TBI","皮层","snRNA")
  rm(hum,hastro); gc(verbose=FALSE) }

# ── 汇总 ──
cat("\n\n══════════ 全数据集星形谷氨酸处理面板 ══════════\n")
tab <- bind_rows(RES)
print(as.data.frame(tab), row.names=FALSE)
write.csv(tab, file.path(FIGOUT,"glutamate_panel_ALL_datasets.csv"), row.names=FALSE)

cat("\n判读：看 GLT-1/GLAST/GS 在轻度→中度→人类 全谱系里是不是一致↓、效应量梯度\n")
cat("\n========= 完成 =========\n")
