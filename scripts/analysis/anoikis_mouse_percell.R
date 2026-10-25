#!/usr/bin/env Rscript
# 补跑: GSE271769 小鼠神经元 anoikis 每细胞数据导出 (供 violin 图)
suppressMessages({library(Seurat); library(Matrix)})
set.seed(1)
RAW <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/raw_data/GSE271769"
OUT <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/anoikis_out"

h5s <- list.files(RAW, pattern="filtered_feature_bc_matrix.h5$", full.names=TRUE)
meta_tp <- function(f){
  b <- basename(f)
  if (grepl("NoTBI",b)) "NoTBI" else if (grepl("_D1_",b)) "D1" else
  if (grepl("_D7_",b)) "D7" else if (grepl("_D28_",b)) "D28" else NA
}
sid <- function(f) sub("_filtered.*","",basename(f))

M <- list(neuron=c("Snap25","Syt1","Rbfox3","Meg3","Syn1"),
          astro=c("Aqp4","Gfap","Slc1a3","Slc1a2","Aldoc","Gja1"),
          micro=c("Cx3cr1","Hexb","C1qa","Ctss","Csf1r"),
          oligo=c("Plp1","Mbp","Mog","Mobp"),
          opc=c("Pdgfra","Cspg4"),vasc=c("Cldn5","Flt1","Pdgfrb","Vtn"))
survival_genes <- c("Ptk2","Akt1","Akt2","Akt3","Bcl2","Mcl1","Bcl2l1","Birc5","Mapk1")
death_genes    <- c("Bcl2l11","Bmf","Bax","Bad","Casp9","Casp3","Cycs","Apaf1")

objs <- list()
for (f in h5s){
  m <- Read10X_h5(f); if (is.list(m)) m <- m[["Gene Expression"]]
  o <- CreateSeuratObject(m, project=sid(f), min.cells=3, min.features=200)
  o$timepoint <- meta_tp(f); o$sample <- sid(f)
  o[["percent.mt"]] <- PercentageFeatureSet(o, pattern="^mt-")
  o <- subset(o, subset = nFeature_RNA>200 & nFeature_RNA<7000 & percent.mt<20)
  o <- NormalizeData(o, verbose=FALSE)
  objs[[sid(f)]] <- o
}
obj <- merge(objs[[1]], objs[-1]); obj <- JoinLayers(obj)

sc <- sapply(names(M), function(k){
  g <- intersect(M[[k]], rownames(obj))
  Matrix::colMeans(GetAssayData(obj, layer="data")[g,,drop=FALSE])
})
obj$celltype <- colnames(sc)[max.col(sc, ties.method="first")]
neu <- subset(obj, subset=celltype=="neuron")
neu <- AddModuleScore(neu, features=list(survival_genes), name="survival_score", search=TRUE)
neu <- AddModuleScore(neu, features=list(death_genes),    name="death_score",    search=TRUE)
neu$survival <- neu$survival_score1
neu$death    <- neu$death_score1
neu$anoikis  <- neu$death - neu$survival

cell_df <- data.frame(sample=neu$sample, timepoint=neu$timepoint,
                      survival=neu$survival, death=neu$death, anoikis=neu$anoikis)
write.csv(cell_df, file.path(OUT,"anoikis_mouse_per_cell.csv"), row.names=FALSE)
cat("小鼠每细胞 n=", nrow(cell_df), "\n")

# 关键基因表达 (heatmap)
key_genes <- c(death_genes, survival_genes)
key_genes <- intersect(key_genes, rownames(neu))
gene_df <- data.frame(t(as.matrix(GetAssayData(neu, layer='data')[key_genes,])))
gene_df$timepoint <- neu$timepoint
gene_df$sample <- neu$sample
write.csv(gene_df, file.path(OUT,"anoikis_mouse_gene_expr.csv"), row.names=FALSE)
cat("DONE -> 小鼠每细胞 + 基因表达已导出\n")
