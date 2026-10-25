#!/usr/bin/env Rscript
# SFig 跨脑区对照 v2 - 直接从 raw h5 读 GSE271769 D7 + NoTBI (避开 v5 RDS OOM)
# 数据: 8 个 h5 (NoTBI rep1-4 + D7 rep1-4), mice 4 vs 4, CCI 重度皮层 7d
set.seed(42)
suppressMessages({library(Seurat); library(Matrix); library(igraph); library(ggplot2)})

RAW <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/raw_data/GSE271769"
OUT <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/upgrade_SFig_CCI_cortex"
PANEL <- c("Nrxn1","Nrxn2","Nrxn3","Negr1","Cntn2","Cntn4","Cntn5",
           "Cadm1","Cadm2","Cadm3","Cdh8","Cdh11","Cdh13","Nlgn1","Nlgn2","Nlgn3")

samples <- list(
  control = c("GSM8384872_NoTBI_rep1","GSM8384873_NoTBI_rep2",
              "GSM8384874_NoTBI_rep3","GSM8384875_NoTBI_rep4"),
  tbi     = c("GSM8384864_D7_rep1","GSM8384865_D7_rep2",
              "GSM8384866_D7_rep3","GSM8384867_D7_rep4"))

cat("[S2.1] read 8 h5 files\n")
obj_list <- list()
for (cond in names(samples)) {
  for (samp in samples[[cond]]) {
    f <- file.path(RAW, paste0(samp, "_filtered_feature_bc_matrix.h5"))
    cat("  read:", samp, "\n")
    m <- Read10X_h5(f)
    s <- CreateSeuratObject(counts=m, project=samp, min.cells=3, min.features=200)
    s$sample <- samp; s$condition <- cond
    obj_list[[samp]] <- s
    rm(m); gc(verbose=FALSE)
  }
}
cat("\n[S2.2] merge\n")
seu <- merge(obj_list[[1]], y=obj_list[-1], add.cell.ids=names(obj_list))
seu <- JoinLayers(seu, assay="RNA")
rm(obj_list); gc(verbose=FALSE)
cat("  total cells after QC:", ncol(seu), "\n")
print(table(seu$condition, seu$sample))

cat("\n[S2.3] cell-level neuron filter (Snap25 > Aqp4/Gfap/Hexb/Plp1)\n")
expr <- GetAssayData(seu, layer="counts", assay="RNA")
mu_cell <- function(g) if (g %in% rownames(expr)) expr[g,] else rep(0, ncol(expr))
sn <- mu_cell("Snap25"); sy <- mu_cell("Syt1")
aq <- mu_cell("Aqp4"); gf <- mu_cell("Gfap"); hx <- mu_cell("Hexb"); pl <- mu_cell("Plp1")
neuron_score <- (sn + sy) / 2
other_score  <- pmax(aq, gf, hx, pl)
is_neuron <- (neuron_score >= 1) & (neuron_score > 2 * other_score)
cat("  neuron-classified:", sum(is_neuron), "/", ncol(seu),
    sprintf(" (%.1f%%)\n", 100*mean(is_neuron)))
neu <- seu[, is_neuron]
rm(seu); gc(verbose=FALSE)
cat("  per-condition neurons:\n"); print(table(neu$condition))
cat("  per-sample neurons:\n"); print(table(neu$sample))

cat("\n[S2.4] gate (entry purity)\n")
expr2 <- GetAssayData(neu, layer="counts", assay="RNA")
mu <- function(g) if (g %in% rownames(expr2)) mean(expr2[g,]) else 0
ratio <- mean(c(mu("Snap25"), mu("Syt1"))) / max(0.01, mean(c(mu("Aqp4"), mu("Gfap"))))
cat(sprintf("  neuron/astro ratio = %.2fx\n", ratio))
if (ratio < 3) cat("  WARN: ratio low (", ratio, ")\n")

cat("\n[S2.5] binarize panel\n")
present <- intersect(PANEL, rownames(expr2))
cat("  panel present:", length(present), "/", length(PANEL), "\n")
mat_bin <- as.matrix(expr2[present, , drop=FALSE] >= 1) * 1
storage.mode(mat_bin) <- "integer"

cat("\n[S2.6] per-group hub network (control vs tbi)\n")
grp <- as.character(neu$condition)
grp_levels <- c("control","tbi")
hub_results <- list(); edge_results <- list(); permR_results <- list()
for (g in grp_levels) {
  cells_g <- which(grp == g)
  cat(sprintf("\n  group %s: %d cells\n", g, length(cells_g)))
  M <- t(mat_bin[, cells_g, drop=FALSE])
  R_obs <- suppressWarnings(cor(M)); R_obs[is.na(R_obs)] <- 0; diag(R_obs) <- 0
  cat("  perm x500 ...\n")
  R_perm <- array(0, dim=c(nrow(R_obs), ncol(R_obs), 500))
  for (k in 1:500) {
    Mp <- apply(M, 2, sample)
    Rp <- suppressWarnings(cor(Mp)); Rp[is.na(Rp)] <- 0; diag(Rp) <- 0
    R_perm[,,k] <- Rp
  }
  q99 <- apply(R_perm, c(1,2), quantile, probs=0.99, na.rm=TRUE); thr <- q99 + 0.05
  adj <- (R_obs > thr) * 1; diag(adj) <- 0
  cat("  edges:", sum(adj)/2, "\n")
  gr <- if (sum(adj)>0) graph_from_adjacency_matrix(adj, mode="undirected") else make_empty_graph(n=ncol(R_obs), directed=FALSE)
  V(gr)$name <- rownames(R_obs)
  memb <- if (sum(adj)>0) membership(cluster_louvain(gr)) else rep(1, ncol(R_obs))
  hub_results[[g]] <- data.frame(gene=rownames(R_obs), hub=as.integer(memb), obs_mean_r=rowMeans(R_obs))
  edges_df <- as.data.frame(as_edgelist(gr))
  if (nrow(edges_df)>0) { colnames(edges_df) <- c("from","to")
    edges_df$r_obs <- mapply(function(a,b) R_obs[a,b], edges_df$from, edges_df$to) }
  edge_results[[g]] <- edges_df
  permR_results[[g]] <- list(R_obs=R_obs, thr=thr)
}

write.csv(do.call(rbind, lapply(grp_levels, function(g) data.frame(group=g, hub_results[[g]]))),
  file.path(OUT, "S_hub_membership.csv"), row.names=FALSE)
summary_df <- data.frame(group=grp_levels,
  n_cells = sapply(grp_levels, function(g) sum(grp==g)),
  n_edges = sapply(grp_levels, function(g) nrow(edge_results[[g]])),
  mean_obsR = sapply(grp_levels, function(g) mean(permR_results[[g]]$R_obs[upper.tri(permR_results[[g]]$R_obs)])))
cat("\n[S2.7] summary\n"); print(summary_df)
write.csv(summary_df, file.path(OUT, "S_summary.csv"), row.names=FALSE)

col2 <- c(control="#4DAECD", tbi="#E07B54")
p1 <- ggplot(summary_df, aes(group, n_edges, fill=group)) + geom_col(width=0.55, colour="grey20") +
  scale_fill_manual(values=col2, guide="none") +
  labs(title="SFig. CCI cortex 7d neuron adhesion co-positivity (GSE271769)",
       subtitle="D7 vs NoTBI, n=4 vs 4 mice — cross-region control",
       y="# edges (r > q99+0.05)") + theme_bw(base_size=11)
ggsave(file.path(OUT,"S_edges_bar.pdf"), p1, width=5, height=4)
ggsave(file.path(OUT,"S_edges_bar.png"), p1, width=5, height=4, dpi=200)

p2 <- ggplot(summary_df, aes(group, mean_obsR, fill=group)) + geom_col(width=0.55, colour="grey20") +
  scale_fill_manual(values=col2, guide="none") +
  labs(title="SFig. Mean co-positivity correlation (CCI cortex 7d)",
       subtitle="GSE271769 D7 vs NoTBI",
       y="mean(R_obs) upper triangle") + theme_bw(base_size=11)
ggsave(file.path(OUT,"S_meanR_bar.pdf"), p2, width=5, height=4)
ggsave(file.path(OUT,"S_meanR_bar.png"), p2, width=5, height=4, dpi=200)

cat("\n[S2.8] DONE at:\n  ", OUT, "\n")
