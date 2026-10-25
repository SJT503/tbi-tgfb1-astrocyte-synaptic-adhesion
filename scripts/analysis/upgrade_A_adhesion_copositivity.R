#!/usr/bin/env Rscript
# upgrade_A_adhesion_copositivity.R — v2.2 升级 A
# Co-positivity hub network on Neuron subset (Sham / CCI_Scr / CCI_KO)
# 范式: SenePy (Sanborn 2025 Nat Commun)
# 输入: all_cells_annotated.rds (自测三组), subset Neuron
# 输出: per-group hub network + hub activity per cell + 三组对比图
set.seed(42)
suppressMessages({
  library(Seurat); library(Matrix); library(igraph); library(ggplot2)
})

# ── PATHS ────────────────────────────────────────────────────────────
RDS  <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/wet_lab_data/results/seurat/all_cells_annotated.rds"
OUT  <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/upgrade_A_copositivity"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
GRP_COL <- "stim"   # values: Sham / CCI_Scr / CCI_TGFbetaR1_KO

# ── ADHESION PANEL (per v3 §4 line 60) ───────────────────────────────
PANEL <- c("Nrxn1","Nrxn2","Nrxn3", "Negr1",
           "Cntn2","Cntn4","Cntn5",
           "Cadm1","Cadm2","Cadm3",
           "Cdh8","Cdh11","Cdh13",
           "Nlgn1","Nlgn2","Nlgn3")  # NCAM 不列入(v2.1 框架 NCAM 不降)

# ── LOAD + GATE (cell-type purity for Neuron subset) ────────────────
cat("[A.1] load + subset Neuron\n")
seu <- readRDS(RDS); DefaultAssay(seu) <- "RNA"
neu <- subset(seu, subset = celltype == "Neuron")
cat("  Neuron cells:", ncol(neu), "\n")
expr <- GetAssayData(neu, layer="counts", assay="RNA")
mu <- function(g) if (g %in% rownames(expr)) mean(expr[g,]) else NA
m  <- c(Snap25=mu("Snap25"), Syt1=mu("Syt1"), Aqp4=mu("Aqp4"), Gfap=mu("Gfap"))
cat("  marker means:", paste(names(m), round(m,3), sep="=", collapse=", "), "\n")
neuron_mean <- mean(m[c("Snap25","Syt1")], na.rm=TRUE)
astro_mean  <- mean(m[c("Aqp4","Gfap")],   na.rm=TRUE)
ratio <- neuron_mean / astro_mean
cat(sprintf("  purity (neuron/astro) = %.2fx (need >=3)\n", ratio))
if (is.na(ratio) || ratio < 3)
  stop(sprintf("GATE FAIL: neuron purity %.2f < 3", ratio))
cat("  [GATE PASS]\n")

# ── BINARIZE PANEL EXPRESSION (>=1 UMI → 1) ─────────────────────────
present <- intersect(PANEL, rownames(expr))
cat("\n[A.2] panel genes present:", length(present), "/", length(PANEL), "\n")
cat("  ", paste(present, collapse=", "), "\n")
if (length(present) < 5) stop("panel coverage too low")
mat_bin <- as.matrix(expr[present, , drop=FALSE] >= 1) * 1
storage.mode(mat_bin) <- "integer"

# ── PER-GROUP CO-POSITIVITY (Pearson + 500 perms) ───────────────────
grp <- as.character(neu@meta.data[[GRP_COL]])
grp_levels <- c("Sham","CCI_Scr","CCI_TGFbetaR1_KO")
hub_results <- list()
edge_results <- list()
permR_results <- list()
for (g in grp_levels) {
  cat(sprintf("\n[A.3] group %s\n", g))
  cells_g <- which(grp == g)
  cat("  cells:", length(cells_g), "\n")
  M <- t(mat_bin[, cells_g, drop=FALSE])   # cells × genes
  # observed pairwise Pearson
  R_obs <- suppressWarnings(cor(M))
  R_obs[is.na(R_obs)] <- 0
  diag(R_obs) <- 0
  # permutation null (shuffle each column independently)
  cat("  permutation x500 ...\n")
  R_perm <- array(0, dim=c(nrow(R_obs), ncol(R_obs), 500))
  for (k in 1:500) {
    Mp <- apply(M, 2, sample)
    Rp <- suppressWarnings(cor(Mp))
    Rp[is.na(Rp)] <- 0; diag(Rp) <- 0
    R_perm[,,k] <- Rp
  }
  q99 <- apply(R_perm, c(1,2), quantile, probs=0.99, na.rm=TRUE)
  thr <- q99 + 0.05
  # adjacency
  adj <- (R_obs > thr) * 1
  diag(adj) <- 0
  cat("  edges (r>q99+0.05):", sum(adj)/2, "\n")
  # louvain hubs
  if (sum(adj) > 0) {
    gr <- graph_from_adjacency_matrix(adj, mode="undirected", weighted=NULL)
    V(gr)$name <- rownames(R_obs)
    comm <- cluster_louvain(gr)
    memb <- membership(comm)
  } else {
    gr <- make_empty_graph(n=ncol(R_obs), directed=FALSE)
    V(gr)$name <- rownames(R_obs)
    memb <- rep(1, ncol(R_obs))
  }
  hub_results[[g]] <- data.frame(gene=rownames(R_obs), hub=as.integer(memb),
                                 obs_mean_r=rowMeans(R_obs))
  edges_df <- as.data.frame(as_edgelist(gr))
  if (nrow(edges_df) > 0) {
    colnames(edges_df) <- c("from","to")
    edges_df$r_obs <- mapply(function(a,b) R_obs[a,b], edges_df$from, edges_df$to)
  }
  edge_results[[g]] <- edges_df
  permR_results[[g]] <- list(R_obs=R_obs, thr=thr)
}

# ── HUB ACTIVITY PER CELL (mean binarized panel within each hub) ─────
cat("\n[A.4] per-cell hub activity (using Sham hub assignment as ref)\n")
ref_hubs <- hub_results[["Sham"]]
hub_act <- list()
for (g in grp_levels) {
  cells_g <- which(grp == g)
  M <- mat_bin[, cells_g, drop=FALSE]
  per_hub <- sapply(sort(unique(ref_hubs$hub)), function(h) {
    genes_h <- ref_hubs$gene[ref_hubs$hub == h]
    genes_h <- intersect(genes_h, rownames(M))
    if (length(genes_h) == 0) return(rep(NA, ncol(M)))
    colMeans(M[genes_h, , drop=FALSE])
  })
  colnames(per_hub) <- paste0("hub", sort(unique(ref_hubs$hub)))
  per_hub <- as.data.frame(per_hub); per_hub$group <- g
  hub_act[[g]] <- per_hub
}
hub_act_all <- do.call(rbind, hub_act)
hub_act_all$group <- factor(hub_act_all$group, levels=grp_levels)

# ── SAVE OUTPUTS ────────────────────────────────────────────────────
write.csv(do.call(rbind, lapply(grp_levels, function(g)
  data.frame(group=g, hub_results[[g]]))),
  file.path(OUT, "A_hub_membership.csv"), row.names=FALSE)
for (g in grp_levels) {
  write.csv(edge_results[[g]],
    file.path(OUT, sprintf("A_edges_%s.csv", g)), row.names=FALSE)
  write.csv(permR_results[[g]]$R_obs,
    file.path(OUT, sprintf("A_R_obs_%s.csv", g)))
}
write.csv(hub_act_all,
  file.path(OUT, "A_hub_activity_percell.csv"), row.names=FALSE)

# ── SUMMARY STATS ───────────────────────────────────────────────────
summary_df <- data.frame(group=grp_levels,
  n_cells = sapply(grp_levels, function(g) sum(grp==g)),
  n_edges = sapply(grp_levels, function(g) nrow(edge_results[[g]])),
  n_hubs  = sapply(grp_levels, function(g) length(unique(hub_results[[g]]$hub))),
  mean_obsR = sapply(grp_levels, function(g) mean(permR_results[[g]]$R_obs[upper.tri(permR_results[[g]]$R_obs)]))
)
print(summary_df)
write.csv(summary_df, file.path(OUT, "A_summary.csv"), row.names=FALSE)

# ── QUICK FIGURES (3-panel) ─────────────────────────────────────────
cat("\n[A.5] quick figures\n")
# (1) edges count bar
p1 <- ggplot(summary_df, aes(group, n_edges, fill=group)) +
  geom_col(width=0.6, colour="grey20") +
  scale_fill_manual(values=c(Sham="#4DAECD", CCI_Scr="#E07B54", CCI_TGFbetaR1_KO="#6A9A5B"),
                    guide="none") +
  labs(title="A. Co-positivity hub edges per group", y="# edges (r>q99+0.05)") +
  theme_bw(base_size=11)
ggsave(file.path(OUT, "A_edges_bar.pdf"), p1, width=5, height=4)
ggsave(file.path(OUT, "A_edges_bar.png"), p1, width=5, height=4, dpi=200)

# (2) hub activity per cell, faceted
hub_long <- reshape(hub_act_all, varying=grep("^hub", colnames(hub_act_all), value=TRUE),
                    v.names="activity", timevar="hub",
                    times=grep("^hub", colnames(hub_act_all), value=TRUE),
                    direction="long")
p2 <- ggplot(hub_long, aes(group, activity, fill=group)) +
  geom_violin(width=0.85, scale="width", trim=TRUE, colour="grey25") +
  facet_wrap(~hub, scales="free_y") +
  scale_fill_manual(values=c(Sham="#4DAECD", CCI_Scr="#E07B54", CCI_TGFbetaR1_KO="#6A9A5B"),
                    guide="none") +
  labs(title="A. Per-cell hub activity by group",
       y="Mean binarized expression in hub") +
  theme_bw(base_size=10) + theme(axis.text.x=element_text(angle=20, hjust=1))
ggsave(file.path(OUT, "A_hub_activity_violins.pdf"), p2, width=7, height=5)
ggsave(file.path(OUT, "A_hub_activity_violins.png"), p2, width=7, height=5, dpi=200)

cat("\n[A.6] DONE. Outputs at:\n  ", OUT, "\n")
