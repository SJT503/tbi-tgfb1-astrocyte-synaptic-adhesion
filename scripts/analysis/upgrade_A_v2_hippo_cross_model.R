#!/usr/bin/env Rscript
# upgrade_A_v2 — Hippocampus neuron hub network across two focal TBI models
# 海马 KO 匹配 + 跨模型保守性 (concussion GSE101901 + FPI GSE180862)
set.seed(42)
suppressMessages({library(Seurat); library(Matrix); library(igraph); library(ggplot2)})

RDS <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/processed_data/integrated/mouse_integrated_annotated.rds"
OUT <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/upgrade_A_v2_hippo"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
PANEL <- c("Nrxn1","Nrxn2","Nrxn3","Negr1","Cntn2","Cntn4","Cntn5",
           "Cadm1","Cadm2","Cadm3","Cdh8","Cdh11","Cdh13","Nlgn1","Nlgn2","Nlgn3")

cat("[A2.1] load + subset hippocampus + neuron\n")
seu <- readRDS(RDS); DefaultAssay(seu) <- "RNA"
hip <- subset(seu, subset = brain_region == "hippocampus" & celltype_auto == "neuron")
cat("  hippocampus neuron cells:", ncol(hip), "\n")
cat("  by condition x dataset x timepoint:\n")
print(xtabs(~ condition + dataset + timepoint, data=hip@meta.data))

cat("\n[A2.2] entry gate (Snap25/Syt1 vs Aqp4/Gfap)\n")
expr <- GetAssayData(hip, layer="counts", assay="RNA")
mu <- function(g) if (g %in% rownames(expr)) mean(expr[g,]) else NA
m <- c(Snap25=mu("Snap25"), Syt1=mu("Syt1"), Aqp4=mu("Aqp4"), Gfap=mu("Gfap"))
cat("  marker means:", paste(names(m), round(m,3), sep="=", collapse=", "), "\n")
ratio <- mean(m[c("Snap25","Syt1")], na.rm=TRUE) / mean(m[c("Aqp4","Gfap")], na.rm=TRUE)
cat(sprintf("  neuron/astro ratio = %.2fx (need >=3)\n", ratio))
if (is.na(ratio) || ratio < 3) stop(sprintf("GATE FAIL: %.2f < 3", ratio))
cat("  [GATE PASS]\n")

cat("\n[A2.3] binarize panel\n")
present <- intersect(PANEL, rownames(expr))
cat("  panel present:", length(present), "/", length(PANEL),"\n")
mat_bin <- as.matrix(expr[present, , drop=FALSE] >= 1) * 1; storage.mode(mat_bin) <- "integer"

cat("\n[A2.4] per-group hub network (control vs tbi; all 3 timepoints combined)\n")
grp <- as.character(hip$condition)
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
  gr <- if (sum(adj)>0) graph_from_adjacency_matrix(adj, mode="undirected", weighted=NULL)
        else make_empty_graph(n=ncol(R_obs), directed=FALSE)
  V(gr)$name <- rownames(R_obs)
  memb <- if (sum(adj)>0) membership(cluster_louvain(gr)) else rep(1, ncol(R_obs))
  hub_results[[g]] <- data.frame(gene=rownames(R_obs), hub=as.integer(memb),
                                 obs_mean_r=rowMeans(R_obs))
  edges_df <- as.data.frame(as_edgelist(gr))
  if (nrow(edges_df)>0) { colnames(edges_df) <- c("from","to")
    edges_df$r_obs <- mapply(function(a,b) R_obs[a,b], edges_df$from, edges_df$to) }
  edge_results[[g]] <- edges_df
  permR_results[[g]] <- list(R_obs=R_obs, thr=thr)
}

# also per-model breakdown (concussion vs FPI within hippocampus)
cat("\n[A2.5] per-model breakdown for cross-model robustness\n")
ds <- as.character(hip$dataset)
model_results <- list()
for (m_name in c("GSE101901","GSE180862")) {
  for (cnd in c("control","tbi")) {
    key <- paste(m_name, cnd, sep="|")
    cells_k <- which(ds == m_name & grp == cnd)
    if (length(cells_k) < 30) { cat("  skip", key, "(", length(cells_k), "cells)\n"); next }
    cat(sprintf("  %s: %d cells\n", key, length(cells_k)))
    M <- t(mat_bin[, cells_k, drop=FALSE])
    R_obs <- suppressWarnings(cor(M)); R_obs[is.na(R_obs)] <- 0; diag(R_obs) <- 0
    model_results[[key]] <- list(R_obs=R_obs, n_cells=length(cells_k),
                                 mean_r=mean(R_obs[upper.tri(R_obs)]))
  }
}

# ── SAVE ────────────────────────────────────────────────────────────
write.csv(do.call(rbind, lapply(grp_levels, function(g) data.frame(group=g, hub_results[[g]]))),
  file.path(OUT, "A2_hub_membership.csv"), row.names=FALSE)
for (g in grp_levels) {
  write.csv(edge_results[[g]], file.path(OUT, sprintf("A2_edges_%s.csv", g)), row.names=FALSE)
  write.csv(permR_results[[g]]$R_obs, file.path(OUT, sprintf("A2_R_obs_%s.csv", g)))
}
summary_df <- data.frame(group=grp_levels,
  n_cells = sapply(grp_levels, function(g) sum(grp==g)),
  n_edges = sapply(grp_levels, function(g) nrow(edge_results[[g]])),
  n_hubs  = sapply(grp_levels, function(g) length(unique(hub_results[[g]]$hub))),
  mean_obsR = sapply(grp_levels, function(g) mean(permR_results[[g]]$R_obs[upper.tri(permR_results[[g]]$R_obs)])))
cat("\n[A2.6] summary\n"); print(summary_df)
write.csv(summary_df, file.path(OUT, "A2_summary.csv"), row.names=FALSE)

# model-level summary
model_df <- data.frame(
  model_cond = names(model_results),
  n_cells = sapply(model_results, function(x) x$n_cells),
  mean_r  = sapply(model_results, function(x) x$mean_r))
write.csv(model_df, file.path(OUT, "A2_per_model_summary.csv"), row.names=FALSE)
cat("\n[A2.7] per-model summary\n"); print(model_df)

# ── FIGURES ─────────────────────────────────────────────────────────
cat("\n[A2.8] figures\n")
col2 <- c(control="#4DAECD", tbi="#E07B54")
# (1) edge count
p1 <- ggplot(summary_df, aes(group, n_edges, fill=group)) +
  geom_col(width=0.55, colour="grey20") +
  scale_fill_manual(values=col2, guide="none") +
  labs(title="A2. Hippocampus neuron adhesion co-positivity edges",
       subtitle="GSE101901 (concussion) + GSE180862 (FPI), all timepoints",
       y="# edges (r > q99 + 0.05)") + theme_bw(base_size=11)
ggsave(file.path(OUT,"A2_edges_bar.pdf"), p1, width=5, height=4)
ggsave(file.path(OUT,"A2_edges_bar.png"), p1, width=5, height=4, dpi=200)

# (2) mean R bar
p2 <- ggplot(summary_df, aes(group, mean_obsR, fill=group)) +
  geom_col(width=0.55, colour="grey20") +
  scale_fill_manual(values=col2, guide="none") +
  labs(title="A2. Mean observed co-positivity correlation",
       y="mean(R_obs) over panel pairs") + theme_bw(base_size=11)
ggsave(file.path(OUT,"A2_meanR_bar.pdf"), p2, width=5, height=4)
ggsave(file.path(OUT,"A2_meanR_bar.png"), p2, width=5, height=4, dpi=200)

# (3) per-model x condition mean R
model_df$model <- gsub("\\|.*","",model_df$model_cond)
model_df$cond  <- gsub(".*\\|","",model_df$model_cond)
p3 <- ggplot(model_df, aes(model, mean_r, fill=cond)) +
  geom_col(position=position_dodge(0.75), width=0.65, colour="grey20") +
  scale_fill_manual(values=col2) +
  labs(title="A2. Per-model robustness (hippocampus)",
       subtitle="Concussion (GSE101901) vs FPI (GSE180862) — does TBI collapse hold cross-model?",
       y="mean R (upper triangle)", x=NULL) + theme_bw(base_size=11)
ggsave(file.path(OUT,"A2_per_model_meanR.pdf"), p3, width=5, height=4)
ggsave(file.path(OUT,"A2_per_model_meanR.png"), p3, width=5, height=4, dpi=200)

cat("\n[A2.9] DONE at:\n  ", OUT, "\n")
