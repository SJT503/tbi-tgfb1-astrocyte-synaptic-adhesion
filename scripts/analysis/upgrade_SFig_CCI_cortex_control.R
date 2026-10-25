#!/usr/bin/env Rscript
# SFig 跨脑区对照: CCI 皮层 7d (GSE271769 D7+NoTBI) hub network
# 验证黏附塌陷 signature 不限于 specific brain region
# 注: v2 severity_integrated 无 celltype 标注,用 Snap25/Syt1 marker 直接细胞级筛 neuron
set.seed(42)
suppressMessages({library(Seurat); library(Matrix); library(igraph); library(ggplot2)})

RDS <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/processed_data/v2/severity_integrated_allcells.rds"
OUT <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/upgrade_SFig_CCI_cortex"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
PANEL <- c("Nrxn1","Nrxn2","Nrxn3","Negr1","Cntn2","Cntn4","Cntn5",
           "Cadm1","Cadm2","Cadm3","Cdh8","Cdh11","Cdh13","Nlgn1","Nlgn2","Nlgn3")

cat("[S.1] load v2 RDS (420k cells) + JoinLayers (Seurat v5)\n")
seu <- readRDS(RDS); DefaultAssay(seu) <- "RNA"
seu <- JoinLayers(seu, assay="RNA")
cat("  layers after join:", paste(Layers(seu, "RNA"), collapse=","), "\n")

cat("\n[S.2] subset GSE271769 D7 + NoTBI (CCI 7d cortex)\n")
sub <- subset(seu, subset = dataset == "GSE271769" & time_point %in% c("D7","NoTBI"))
cat("  cells:", ncol(sub), "\n")
cat("  by time_point:\n"); print(table(sub$time_point))
rm(seu); gc(verbose=FALSE)

cat("\n[S.3] marker-based neuron filter (Snap25 + Syt1 cell-level)\n")
expr <- GetAssayData(sub, layer="counts", assay="RNA")
m_present <- intersect(c("Snap25","Syt1","Aqp4","Gfap","Hexb","Csf1r","Plp1","Mbp"), rownames(expr))
cat("  markers present:", paste(m_present, collapse=","), "\n")
sn <- if ("Snap25" %in% rownames(expr)) expr["Snap25",] else rep(0, ncol(expr))
sy <- if ("Syt1"   %in% rownames(expr)) expr["Syt1",]   else rep(0, ncol(expr))
aq <- if ("Aqp4"   %in% rownames(expr)) expr["Aqp4",]   else rep(0, ncol(expr))
gf <- if ("Gfap"   %in% rownames(expr)) expr["Gfap",]   else rep(0, ncol(expr))
hx <- if ("Hexb"   %in% rownames(expr)) expr["Hexb",]   else rep(0, ncol(expr))
pl <- if ("Plp1"   %in% rownames(expr)) expr["Plp1",]   else rep(0, ncol(expr))
neuron_score <- (sn + sy) / 2
other_score  <- pmax(aq, gf, hx, pl)
is_neuron <- (neuron_score >= 1) & (neuron_score > 2 * other_score)
cat("  neuron-classified cells:", sum(is_neuron), "/", ncol(sub),
    sprintf("(%.1f%%)\n", 100*mean(is_neuron)))
neu <- sub[, is_neuron]
cat("  per timepoint:\n"); print(table(neu$time_point))

cat("\n[S.4] gate (entry purity)\n")
expr2 <- GetAssayData(neu, layer="counts", assay="RNA")
mu <- function(g) if (g %in% rownames(expr2)) mean(expr2[g,]) else 0
ratio <- mean(c(mu("Snap25"), mu("Syt1"))) / max(0.01, mean(c(mu("Aqp4"), mu("Gfap"))))
cat(sprintf("  neuron/astro ratio = %.2fx\n", ratio))
if (ratio < 3) { cat("  WARN: ratio low, neuron filter may be loose\n") }

cat("\n[S.5] binarize panel + hub network per condition\n")
present <- intersect(PANEL, rownames(expr2))
cat("  panel present:", length(present), "/", length(PANEL), "\n")
mat_bin <- as.matrix(expr2[present, , drop=FALSE] >= 1) * 1
storage.mode(mat_bin) <- "integer"

# group: D7 = TBI, NoTBI = control
grp <- ifelse(neu$time_point == "D7", "tbi", "control")
grp_levels <- c("control","tbi")
hub_results <- list(); edge_results <- list(); permR_results <- list()
for (g in grp_levels) {
  cells_g <- which(grp == g)
  cat(sprintf("\n  group %s: %d cells\n", g, length(cells_g)))
  if (length(cells_g) < 30) { cat("  skip (too few)\n"); next }
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

write.csv(do.call(rbind, lapply(grp_levels, function(g) {
  if (is.null(hub_results[[g]])) NULL else data.frame(group=g, hub_results[[g]])})),
  file.path(OUT, "S_hub_membership.csv"), row.names=FALSE)
summary_df <- data.frame(group=grp_levels,
  n_cells = sapply(grp_levels, function(g) sum(grp==g, na.rm=TRUE)),
  n_edges = sapply(grp_levels, function(g) if (is.null(edge_results[[g]])) NA else nrow(edge_results[[g]])),
  mean_obsR = sapply(grp_levels, function(g) if (is.null(permR_results[[g]])) NA else mean(permR_results[[g]]$R_obs[upper.tri(permR_results[[g]]$R_obs)])))
cat("\n[S.6] summary\n"); print(summary_df)
write.csv(summary_df, file.path(OUT, "S_summary.csv"), row.names=FALSE)

# figures
col2 <- c(control="#4DAECD", tbi="#E07B54")
p1 <- ggplot(summary_df, aes(group, n_edges, fill=group)) + geom_col(width=0.55, colour="grey20") +
  scale_fill_manual(values=col2, guide="none") +
  labs(title="SFig. CCI cortex 7d neuron adhesion co-positivity edges",
       subtitle="GSE271769 D7 vs NoTBI cortex — cross-region control",
       y="# edges (r > q99+0.05)") + theme_bw(base_size=11)
ggsave(file.path(OUT,"S_edges_bar.pdf"), p1, width=5, height=4)
ggsave(file.path(OUT,"S_edges_bar.png"), p1, width=5, height=4, dpi=200)

p2 <- ggplot(summary_df, aes(group, mean_obsR, fill=group)) + geom_col(width=0.55, colour="grey20") +
  scale_fill_manual(values=col2, guide="none") +
  labs(title="SFig. Mean co-positivity correlation (CCI cortex 7d)",
       y="mean(R_obs)") + theme_bw(base_size=11)
ggsave(file.path(OUT,"S_meanR_bar.pdf"), p2, width=5, height=4)
ggsave(file.path(OUT,"S_meanR_bar.png"), p2, width=5, height=4, dpi=200)

cat("\n[S.7] DONE at:\n  ", OUT, "\n")
