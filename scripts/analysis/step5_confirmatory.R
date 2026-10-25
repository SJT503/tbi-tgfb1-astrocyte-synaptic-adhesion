#!/usr/bin/env Rscript
# ============================================================
# Step 5 — Confirmatory analysis for the Fig 6F mechanism model
#   5a. TGFBR1 / TGFB1 distribution by cell type (mouse + human)
#   5b. astrocyte->neuron cell-cell communication: death vs synaptic axes
#   5c. spatial co-localization: astrocyte TGF-b signalling vs neuronal death
# Data: public datasets only (real processed objects, no simulation)
# ============================================================
suppressMessages({library(Seurat); library(Matrix)})
PD  <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/processed_data"
OUT <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/step5_out"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

expr_summary <- function(seu, genes, group_col, celltype_col, assay="RNA"){
  genes <- genes[genes %in% rownames(seu)]
  m <- GetAssayData(seu, assay=assay, layer="data")
  ct <- as.character(seu@meta.data[[celltype_col]])
  gr <- as.character(seu@meta.data[[group_col]])
  res <- list()
  for(g in genes) for(c in unique(ct)) for(p in unique(gr)){
    idx <- which(ct==c & gr==p)
    if(length(idx)<10) next
    v <- m[g, idx]
    res[[length(res)+1]] <- data.frame(gene=g, celltype=c, group=p,
      n=length(idx), mean=round(mean(v),4), det_pct=round(100*mean(v>0),2))
  }
  do.call(rbind, res)
}

# ---------- PURITY GATE ----------
purity_gate <- function(seu, celltype_col, assay="RNA"){
  m <- GetAssayData(seu, assay=assay, layer="data")
  ct <- as.character(seu@meta.data[[celltype_col]])
  astro <- grep("astro", unique(ct), ignore.case=TRUE, value=TRUE)
  neur  <- grep("neuron", unique(ct), ignore.case=TRUE, value=TRUE)
  if(length(astro) && all(c("Aqp4","Snap25") %in% rownames(m))){
    ai <- which(ct %in% astro)
    r <- mean(m["Aqp4",ai]) / (mean(m["Snap25",ai])+1e-6)
    cat(sprintf("  PURITY astrocyte Aqp4/Snap25 = %.2f %s\n", r, ifelse(r>=3,"PASS","** CHECK **")))
  }
  if(length(neur) && all(c("Aqp4","Snap25") %in% rownames(m))){
    ni <- which(ct %in% neur)
    r <- mean(m["Snap25",ni]) / (mean(m["Aqp4",ni])+1e-6)
    cat(sprintf("  PURITY neuron Snap25/Aqp4 = %.2f %s\n", r, ifelse(r>=3,"PASS","** CHECK **")))
  }
}

cat("\n##### 5a — TGFBR1/TGFB1 by cell type #####\n")
cat("--- MOUSE integrated atlas (GSE101901+GSE180862) ---\n")
mo <- readRDS(file.path(PD,"integrated/mouse_integrated_annotated.rds"))
DefaultAssay(mo) <- "RNA"
purity_gate(mo, "celltype_auto")
mg <- c("Tgfb1","Tgfb2","Tgfb3","Tgfbr1","Tgfbr2")
mres <- expr_summary(mo, mg, "condition", "celltype_auto")
write.csv(mres, file.path(OUT,"5a_mouse_tgfb_byCelltype.csv"), row.names=FALSE)
cat("TGFBR1 (Tgfbr1) by cell type:\n")
print(mres[mres$gene=="Tgfbr1", c("celltype","group","n","mean","det_pct")], row.names=FALSE)
cat("\nTGFB1 (Tgfb1) by cell type:\n")
print(mres[mres$gene=="Tgfb1", c("celltype","group","n","mean","det_pct")], row.names=FALSE)

cat("\n--- HUMAN GSE209552 (pre-computed, scripts/111) ---\n")
hcsv <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/figures/human_tgfb_byCelltype.csv"
if(file.exists(hcsv)){
  h <- read.csv(hcsv)
  cat("TGFBR1 by cell type (human):\n")
  print(h[h$gene=="TGFBR1", c("celltype","group","mean","det")], row.names=FALSE)
} else cat("  human CSV not found\n")

cat("\n##### 5b — astrocyte->neuron communication #####\n")
cc <- readRDS(file.path(PD,"communication/cellchat_like_comparison.rds"))
an <- cc[cc$source=="astrocyte" & cc$target=="neuron",]
# aggregate by pathway
agg <- aggregate(cbind(score_control,score_tbi,delta_tbi_minus_control)~pathway_name, data=an, FUN=sum)
agg$log2FC <- round(log2((agg$score_tbi+1)/(agg$score_control+1)),3)
agg <- agg[order(-agg$delta_tbi_minus_control),]
# classify
synaptic <- c("NRXN","NCAM","NEGR","NGL","CADM","CNTN","CDH","NLG","JAM","NRG","SLITRK")
deathcyto <- c("TNF","TWEAK","FASLG","TRAIL","IL1","IL6","CCL","CXCL","CSF","COMPLEMENT",
               "GRN","SPP1","IL2","IFN","TGFb","LT","GAS")
trophic <- c("NGF","BDNF","GDNF","IGF","FGF","PDGF","VEGF","EGF","PTN","MK","HBEGF","WNT","NT")
cls <- function(p){
  if(p %in% synaptic) return("synaptic/adhesion")
  if(any(sapply(deathcyto, function(x) grepl(x,p,ignore.case=TRUE)))) return("death/cytokine")
  if(p %in% trophic) return("trophic/growth")
  "other/matrix"
}
agg$category <- sapply(agg$pathway_name, cls)
write.csv(agg, file.path(OUT,"5b_astro_to_neuron_pathways.csv"), row.names=FALSE)
cat("ALL astrocyte->neuron pathways (control vs TBI, sorted by delta):\n")
print(agg[,c("pathway_name","category","score_control","score_tbi","log2FC")], row.names=FALSE)
cat("\nTop UP-in-TBI by category:\n")
for(ct in unique(agg$category)){
  s <- agg[agg$category==ct,]
  cat(sprintf("  [%s] n_pathways=%d  total_delta=%.3g\n", ct, nrow(s), sum(s$delta_tbi_minus_control)))
}

cat("\n##### 5c — spatial co-localization (GSE226208 3dpi) #####\n")
sp <- readRDS(file.path(PD,"spatial/visium_tbi_3dpi.rds"))
DefaultAssay(sp) <- "Spatial"
sp <- NormalizeData(sp, verbose=FALSE)
sm <- GetAssayData(sp, assay="Spatial", layer="data")
hav <- function(gs) gs[gs %in% rownames(sm)]
tgfb_sig  <- hav(c("Tgfb1","Tgfbr1","Tgfbr2","Smad2","Smad3","Serpine1","Smad7"))
astro_sig <- hav(c("Aqp4","Gfap","Slc1a3","Aldh1l1","Vim"))
apop_sig  <- hav(c("Bax","Casp3","Casp8","Casp9","Apaf1","Bbc3","Pmaip1"))
necro_sig <- hav(c("Ripk1","Ripk3","Mlkl","Zbp1"))
neuron_sig<- hav(c("Snap25","Rbfox3","Syt1","Meg3"))
score1 <- function(gs) if(length(gs)) Matrix::colMeans(sm[gs,,drop=FALSE]) else rep(NA,ncol(sm))
S <- data.frame(
  tgfb  = score1(tgfb_sig),
  astro = score1(astro_sig),
  apop  = score1(apop_sig),
  necro = score1(necro_sig),
  neuron= score1(neuron_sig))
write.csv(S, file.path(OUT,"5c_spatial_spot_scores.csv"), row.names=TRUE)
cat("Genes used:\n")
cat("  TGF-b signalling:", paste(tgfb_sig,collapse=","), "\n")
cat("  astrocyte:", paste(astro_sig,collapse=","), "\n")
cat("  apoptosis:", paste(apop_sig,collapse=","), "\n")
cat("  necroptosis:", paste(necro_sig,collapse=","), "\n")
cm <- round(cor(S, use="pairwise.complete.obs", method="spearman"),3)
cat("\nSpearman correlation across", ncol(sm), "TBI spots:\n")
print(cm)
# astro-state score columns already in object?
mdc <- colnames(sp@meta.data)
stcols <- grep("score_[0-9]$|signature", mdc, value=TRUE)
if(length(stcols)){
  cat("\nPre-computed astro-state scores correlation with death scores:\n")
  S2 <- cbind(S[,c("apop","necro")], sp@meta.data[, stcols, drop=FALSE])
  print(round(cor(S2, use="pairwise.complete.obs", method="spearman"),3))
}
cat("\nDONE step5\n")
