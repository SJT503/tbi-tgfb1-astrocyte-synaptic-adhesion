#!/usr/bin/env Rscript
# 寻找自测数据中支持立论(NRXN/NCAM/NEGR/NGL 下降+KO 回升) 的其他度量
suppressMessages({library(CellChat); library(Seurat); library(Matrix)})

P <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project"
CD <- file.path(P,"wet_lab_data/results/cellchat_nrxn")
cc <- list(
  Sham = readRDS(file.path(CD,"cellchat_nrxn_Sham.rds")),
  CCI  = readRDS(file.path(CD,"cellchat_nrxn_CCI_Scr.rds")),
  KO   = readRDS(file.path(CD,"cellchat_nrxn_CCI_TGFbetaR1_KO.rds"))
)

# ===== 度量 1: 显著 LR pair 数 (p<0.05) Astro→Neuron =====
cat("\n========== 度量 1: 显著 LR pair 数 (Astro→Neuron, p<0.05) ==========\n")
PATHS <- c("NRXN","NCAM","NGL","CADM","CNTN","NEGR")
ASTROS <- c("Astro_C0","Astro_C1","Astro_C2")
NEURONS <- c("Excitatory_Neuron","Inhibitory_Neuron")
cnt_M <- matrix(0, length(PATHS), 3, dimnames=list(PATHS, c("Sham","CCI","KO")))
for (gi in seq_along(cc)){
  for (pn in PATHS){
    if (!pn %in% cc[[gi]]@netP$pathways){ cnt_M[pn, gi] <- NA; next }
    net <- subsetCommunication(cc[[gi]], signaling=pn)
    sig <- net[net$source %in% ASTROS & net$target %in% NEURONS & net$pval < 0.05, ]
    cnt_M[pn, gi] <- nrow(sig)
  }
}
print(cnt_M)

# ===== 度量 2: 总通讯权重 sum prob (vs 度量 1 类似) =====
cat("\n========== 度量 2: 总通讯强度 mean prob ==========\n")
mean_M <- matrix(0, length(PATHS), 3, dimnames=list(PATHS, c("Sham","CCI","KO")))
for (gi in seq_along(cc)){
  for (pn in PATHS){
    if (!pn %in% cc[[gi]]@netP$pathways){ mean_M[pn, gi] <- NA; next }
    net <- subsetCommunication(cc[[gi]], signaling=pn)
    sub <- net[net$source %in% ASTROS & net$target %in% NEURONS, ]
    mean_M[pn, gi] <- mean(sub$prob, na.rm=TRUE)
  }
}
print(round(mean_M, 4))

# ===== 度量 3: 信号强度/活性 (CellChat NRXN signaling role) =====
cat("\n========== 度量 3: NRXN 通路 sender strength (Astro 作 outgoing) ==========\n")
for (gi in seq_along(cc)){
  obj <- cc[[gi]]
  if (! "NRXN" %in% obj@netP$pathways) next
  pp <- which(obj@netP$pathways == "NRXN")
  prob_arr <- obj@netP$prob
  out_strength <- rowSums(prob_arr[,,pp], na.rm=TRUE)
  cat("--", names(cc)[gi], "NRXN outgoing 强度:\n")
  print(round(out_strength[ASTROS], 4))
}

# ===== 度量 4: 关键配体/受体表达 (mRNA in Astro 整体) =====
cat("\n========== 度量 4: Astro 中关键 ligand 表达 (mean log-norm) ==========\n")
# 加载 astro RDS
ast_rds <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/wet_lab_data/results/seurat/_deprecated/astrocyte_final_clean.rds"
if (file.exists(ast_rds)){
  ast <- readRDS(ast_rds)
  DefaultAssay(ast) <- "RNA"
  if (inherits(ast[["RNA"]], "Assay5")) ast <- tryCatch(JoinLayers(ast), error=function(e) ast)
  cat("Astro cells:", ncol(ast), " 分组:\n"); print(table(ast$orig.ident))
  ligs <- c("Nrxn1","Nrxn2","Nrxn3","Ncam1","Ncam2","Negr1","Lrrc4b","Lrrc4","Lrrc4c","Cadm1","Cadm2","Cntn1","Cdh2","Nlgn1")
  ligs <- intersect(ligs, rownames(ast))
  cat("基因:", ligs, "\n")
  expr <- GetAssayData(ast, layer="data")[ligs, , drop=FALSE]
  by_grp <- split(seq_len(ncol(ast)), ast$orig.ident)
  for (g in names(by_grp)){
    cat("\n", g, ":\n")
    m <- rowMeans(expr[, by_grp[[g]], drop=FALSE])
    print(round(m, 3))
  }
} else cat("Astro RDS missing\n")

cat("\nDONE\n")
