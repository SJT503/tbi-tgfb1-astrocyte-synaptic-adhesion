#!/usr/bin/env Rscript
# Candidate 2b — human astrocyte->neuron adhesion CellChat, STRICT annotation.
# Fix vs cand2: cluster-level argmax over-called astrocytes (purity 1.8). Here only
# astrocyte-argmax clusters whose AQP4/SNAP25 ratio >= 3 count as Astrocyte.
suppressMessages({library(Seurat); library(CellChat); library(dplyr)})
options(future.globals.maxSize=4e9)
BASE <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project"
HUM  <- file.path(BASE,"wet_lab_data/results/human_validation_v2/processed.rds")
OUT  <- file.path(BASE,"JNI/analysis/cand_out")

MK <- list(Neuron=c("SNAP25","SYT1","RBFOX3"), Astrocyte=c("AQP4","GFAP","SLC1A2"),
  Oligo=c("PLP1","MBP","MOBP"), Microglia=c("CSF1R","P2RY12","CX3CR1"),
  OPC=c("PDGFRA","CSPG4"), Endo=c("CLDN5","FLT1"))

hum <- readRDS(HUM); DefaultAssay(hum) <- "RNA"
hum$condition <- ifelse(hum$condition=="CTL","control","tbi")
ex <- GetAssayData(hum, layer="data", assay="RNA")
cls <- sort(unique(as.character(hum$seurat_clusters)))

# per-cluster: argmax celltype + astrocyte-purity ratio
score <- sapply(names(MK), function(ct){ g<-intersect(MK[[ct]],rownames(ex))
  sapply(cls, function(cl) mean(as.matrix(ex[g, colnames(hum)[as.character(hum$seurat_clusters)==cl],drop=FALSE]))) })
argmax <- colnames(score)[apply(score,1,which.max)]
aqp4 <- ex["AQP4",]; snap <- ex["SNAP25",]
ratio <- sapply(cls, function(cl){ i<-as.character(hum$seurat_clusters)==cl
  mean(aqp4[i])/(mean(snap[i])+1e-6) })
diag <- data.frame(cluster=cls, argmax=argmax, n=as.integer(table(hum$seurat_clusters)[cls]),
                   AQP4_SNAP25_ratio=round(ratio,2))
cat("=== per-cluster annotation + astrocyte purity ===\n"); print(diag, row.names=FALSE)

# STRICT: astrocyte = argmax-astro clusters with ratio>=3; impure astro-argmax -> Mixed(drop)
celltype <- setNames(argmax, cls)
impure <- cls[argmax=="Astrocyte" & ratio<3]
if(length(impure)){ cat("\n⚠ impure astrocyte-argmax clusters dropped:", paste(impure,collapse=","), "\n")
  celltype[impure] <- "Mixed" }
hum$celltype <- unname(celltype[as.character(hum$seurat_clusters)])
hum <- hum[, hum$celltype != "Mixed"]
cat("\ncelltype x condition (strict):\n"); print(table(hum$celltype, hum$condition))
m <- GetAssayData(hum,layer="data",assay="RNA")
ai<-which(hum$celltype=="Astrocyte"); ni<-which(hum$celltype=="Neuron")
cat(sprintf("PURITY astrocyte AQP4/SNAP25=%.1f %s | neuron SNAP25/AQP4=%.1f\n",
  mean(m["AQP4",ai])/(mean(m["SNAP25",ai])+1e-6),
  ifelse(mean(m["AQP4",ai])/(mean(m["SNAP25",ai])+1e-6)>=3,"PASS","**still low**"),
  mean(m["SNAP25",ni])/(mean(m["AQP4",ni])+1e-6)))

runCC <- function(seu){
  dat <- GetAssayData(seu, assay="RNA", layer="data")
  meta <- data.frame(celltype=as.character(seu$celltype), row.names=colnames(seu))
  cc <- createCellChat(object=dat, meta=meta, group.by="celltype")
  cc@DB <- CellChatDB.human
  cc <- subsetData(cc); cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc); cc <- computeCommunProb(cc)
  cc <- filterCommunication(cc, min.cells=10); cc <- computeCommunProbPathway(cc)
  aggregateNet(cc)
}
a2n <- function(cc,tag){ p<-cc@netP$prob
  if(!("Astrocyte"%in%dimnames(p)[[1]])||!("Neuron"%in%dimnames(p)[[2]])) return(NULL)
  v<-p["Astrocyte","Neuron",]; data.frame(pathway=names(v),prob=as.numeric(v),cond=tag) }

res <- list()
for(cd in c("control","tbi")){
  sub <- subset(hum, subset=condition==cd)
  cat("\n--- CellChat:",cd," cells:",ncol(sub),"---\n"); flush.console()
  cc <- tryCatch(runCC(sub), error=function(e){cat("FAIL:",conditionMessage(e),"\n");NULL})
  if(!is.null(cc)){ saveRDS(cc,file.path(OUT,paste0("cand2b_cellchat_",cd,".rds"))); res[[cd]]<-a2n(cc,cd) }
}
tab <- bind_rows(res)
if(nrow(tab)>0){
  w <- reshape(tab, idvar="pathway", timevar="cond", direction="wide")
  colnames(w) <- sub("prob\\.","",colnames(w))
  for(c0 in c("control","tbi")) if(!c0%in%colnames(w)) w[[c0]]<-0
  w[is.na(w)]<-0; w$delta<-round(w$tbi-w$control,6)
  w$pct_change <- round(ifelse(w$control>0,(w$tbi-w$control)/w$control*100,NA),1)
  adh <- c("NRXN","NCAM","NEGR","NGL","CNTN","CADM","NRG","CDH")
  w$adhesion <- w$pathway %in% adh
  w <- w[order(-w$adhesion,-w$control),]
  write.csv(w, file.path(OUT,"cand2b_human_astro2neuron_pathways.csv"), row.names=FALSE)
  cat("\n=== HUMAN (strict) astrocyte->neuron ADHESION pathways ===\n")
  print(w[w$adhesion,c("pathway","control","tbi","pct_change")], row.names=FALSE)
}
cat("\nDONE cand2b\n")
