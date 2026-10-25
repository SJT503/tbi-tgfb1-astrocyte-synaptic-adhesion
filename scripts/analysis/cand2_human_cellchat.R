#!/usr/bin/env Rscript
# Candidate 2 — human astrocyte->neuron adhesion communication (GSE209552),
# TBI vs control. Cross-species confirmation of the 6a adhesion-axis collapse.
suppressMessages({library(Seurat); library(CellChat); library(dplyr)})
options(future.globals.maxSize=4e9)
BASE <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project"
HUM  <- file.path(BASE,"wet_lab_data/results/human_validation_v2/processed.rds")
OUT  <- file.path(BASE,"JNI/analysis/cand_out")
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

MK <- list(Neuron=c("SNAP25","SYT1","RBFOX3"), Astrocyte=c("AQP4","GFAP","SLC1A2"),
  Oligo=c("PLP1","MBP","MOBP"), Microglia=c("CSF1R","P2RY12","CX3CR1"),
  OPC=c("PDGFRA","CSPG4"), Endo=c("CLDN5","FLT1"))

hum <- readRDS(HUM)
DefaultAssay(hum) <- "RNA"
hum$condition <- ifelse(grepl("HuBrainCTL",hum$sample_id),"control",
                 ifelse(grepl("MJ_TBI",hum$sample_id),"tbi",NA))
hum <- hum[, !is.na(hum$condition)]
ex <- GetAssayData(hum, layer="data", assay="RNA")
cls <- sort(unique(hum$seurat_clusters))
sc <- sapply(names(MK), function(ct){ g <- intersect(MK[[ct]], rownames(ex))
  sapply(cls, function(cl) mean(as.matrix(ex[g, colnames(hum)[hum$seurat_clusters==cl],drop=FALSE]))) })
cl2ct <- setNames(colnames(sc)[apply(sc,1,which.max)], cls)
hum$celltype <- unname(cl2ct[as.character(hum$seurat_clusters)])
cat("human cells:", ncol(hum), "\n")
cat("celltype x condition:\n"); print(table(hum$celltype, hum$condition))
# purity gate
pg <- function(s){ m<-GetAssayData(s,layer="data",assay="RNA")
  ai<-which(s$celltype=="Astrocyte"); ni<-which(s$celltype=="Neuron")
  if(length(ai)>10) cat(sprintf("  astro AQP4/SNAP25=%.1f\n", mean(m["AQP4",ai])/(mean(m["SNAP25",ai])+1e-6)))
  if(length(ni)>10) cat(sprintf("  neuron SNAP25/AQP4=%.1f\n", mean(m["SNAP25",ni])/(mean(m["AQP4",ni])+1e-6))) }
pg(hum)

runCC <- function(seu){
  dat <- tryCatch(GetAssayData(seu, assay="RNA", layer="data"),
                  error=function(e) GetAssayData(seu, assay="RNA", slot="data"))
  meta <- data.frame(celltype=as.character(seu$celltype),
                     row.names=colnames(seu), stringsAsFactors=FALSE)
  cc <- createCellChat(object=dat, meta=meta, group.by="celltype")
  cc@DB <- CellChatDB.human
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc)
  cc <- filterCommunication(cc, min.cells=10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc
}

a2n <- function(cc, tag){
  p <- cc@netP$prob               # [source, target, pathway]
  if(!("Astrocyte" %in% dimnames(p)[[1]]) || !("Neuron" %in% dimnames(p)[[2]])){
    cat(sprintf("  [%s] Astrocyte/Neuron missing in netP\n", tag)); return(NULL) }
  v <- p["Astrocyte","Neuron",]
  data.frame(pathway=names(v), prob=as.numeric(v), cond=tag)
}

res <- list()
for(cd in c("control","tbi")){
  sub <- subset(hum, subset = condition == cd)
  cat("\n--- CellChat:", cd, " cells:", ncol(sub), "---\n"); flush.console()
  cc <- tryCatch(runCC(sub), error=function(e){cat("FAIL:",conditionMessage(e),"\n");NULL})
  if(!is.null(cc)){ saveRDS(cc, file.path(OUT, paste0("cand2_cellchat_",cd,".rds")))
    res[[cd]] <- a2n(cc, cd) }
}
tab <- bind_rows(res)
if(nrow(tab)>0){
  w <- reshape(tab, idvar="pathway", timevar="cond", direction="wide")
  colnames(w) <- sub("prob\\.","",colnames(w))
  for(c0 in c("control","tbi")) if(!c0 %in% colnames(w)) w[[c0]] <- 0
  w[is.na(w)] <- 0
  w$delta <- round(w$tbi - w$control, 6)
  w$log2FC <- round(log2((w$tbi+1e-9)/(w$control+1e-9)), 3)
  adh <- c("NRXN","NCAM","NEGR","NGL","CNTN","CADM","NRG","CDH")
  w$adhesion <- w$pathway %in% adh
  w <- w[order(-w$adhesion, w$pathway),]
  write.csv(w, file.path(OUT,"cand2_human_astro2neuron_pathways.csv"), row.names=FALSE)
  cat("\n=== HUMAN astrocyte->neuron pathways (control vs TBI) ===\n")
  cat("--- adhesion pathways ---\n")
  print(w[w$adhesion, c("pathway","control","tbi","log2FC")], row.names=FALSE)
}
cat("\nDONE cand2\n")
