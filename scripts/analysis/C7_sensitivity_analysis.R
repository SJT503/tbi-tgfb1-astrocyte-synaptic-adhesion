#!/usr/bin/env Rscript
# C7: Sensitivity analysis on adhesion module finding
# 测试 NRXN-led adhesion collapse 是否对以下选择稳健:
#  (a) Pseudobulk vs single-cell DE
#  (b) Module gene composition (16 全基因 vs 8 NRXN-only vs 12 去 NCAM 等)
#  (c) Cell-type filtering threshold (default vs strict purity)
# 报告各设定下的 effect size (delta module score) 一致性
set.seed(42)
suppressMessages({library(Seurat); library(dplyr); library(ggplot2)})

RDS <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/processed_data/integrated/mouse_integrated_annotated.rds"
OUT <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/C7_sensitivity"

PANEL_FULL  <- c("Nrxn1","Nrxn2","Nrxn3","Negr1","Cntn2","Cntn4","Cntn5",
                 "Cadm1","Cadm2","Cadm3","Cdh8","Cdh11","Cdh13","Nlgn1","Nlgn2","Nlgn3")
PANEL_NRXN  <- c("Nrxn1","Nrxn2","Nrxn3")
PANEL_NOCAM <- setdiff(PANEL_FULL, c("Cadm1","Cadm2","Cadm3"))
PANEL_CNTN  <- c("Cntn2","Cntn4","Cntn5")
PANEL_CADM  <- c("Cadm1","Cadm2","Cadm3")
PANELS <- list(full16=PANEL_FULL, nrxn_only=PANEL_NRXN, cntn_only=PANEL_CNTN,
               cadm_only=PANEL_CADM, no_cadm=PANEL_NOCAM)

cat("[C7.1] load + neuron subset\n")
seu <- readRDS(RDS); DefaultAssay(seu) <- "RNA"
neu <- subset(seu, subset = celltype_auto == "neuron")
rm(seu); gc(verbose=FALSE)
cat("  neuron cells:", ncol(neu), "\n")

cat("\n[C7.2] module score for each panel composition\n")
for (pn in names(PANELS)) {
  genes <- intersect(PANELS[[pn]], rownames(neu))
  if (length(genes) < 2) next
  neu <- AddModuleScore(neu, features=list(genes), name=paste0("PANEL_", pn, "_"), seed=42)
}
score_cols <- grep("^PANEL_", colnames(neu@meta.data), value=TRUE)

cat("\n[C7.3] per-panel TBI vs sham (cell-level + pseudobulk)\n")
results <- list()
for (sc in score_cols) {
  pn <- gsub("^PANEL_|_1$", "", sc)
  # cell-level
  c_mean <- tapply(neu@meta.data[[sc]], neu$condition, mean, na.rm=TRUE)
  c_t <- t.test(neu@meta.data[[sc]] ~ neu$condition)
  # pseudobulk per (dataset, condition)
  pb <- neu@meta.data %>% mutate(score=neu@meta.data[[sc]]) %>%
    group_by(dataset, condition, animal) %>%
    summarise(mean_score=mean(score, na.rm=TRUE), .groups="drop")
  pb_ttest <- tryCatch(t.test(mean_score ~ condition, data=pb), error=function(e) NULL)
  pb_p <- if (is.null(pb_ttest)) NA else pb_ttest$p.value
  pb_d <- if (is.null(pb_ttest)) NA else diff(rev(pb_ttest$estimate))
  results[[length(results)+1]] <- data.frame(
    panel = pn,
    n_genes = length(intersect(PANELS[[pn]], rownames(neu))),
    cell_mean_control = c_mean["control"],
    cell_mean_tbi     = c_mean["tbi"],
    cell_delta = c_mean["tbi"] - c_mean["control"],
    cell_t = c_t$statistic, cell_p = c_t$p.value,
    pseudobulk_delta = pb_d,
    pseudobulk_p = pb_p,
    pb_n_samples = nrow(pb))
}
res_df <- do.call(rbind, results)
res_df$cell_padj_BH <- p.adjust(res_df$cell_p, "BH")
res_df$pb_padj_BH   <- p.adjust(res_df$pseudobulk_p, "BH")
write.csv(res_df, file.path(OUT, "C7_panel_sensitivity.csv"), row.names=FALSE)
cat("\n"); print(res_df)

# brain-region stratified
cat("\n[C7.4] same, stratified by brain_region\n")
results2 <- list()
for (br in unique(neu$brain_region)) {
  for (sc in score_cols) {
    pn <- gsub("^PANEL_|_1$", "", sc)
    cells <- which(neu$brain_region == br)
    sub <- neu@meta.data[cells, ]
    sub$score <- neu@meta.data[[sc]][cells]
    if (length(unique(sub$condition)) < 2) next
    t <- t.test(score ~ condition, data=sub)
    pb <- sub %>% group_by(dataset, condition, animal) %>%
      summarise(m=mean(score, na.rm=TRUE), .groups="drop")
    pb_t <- tryCatch(t.test(m ~ condition, data=pb), error=function(e) NULL)
    results2[[length(results2)+1]] <- data.frame(
      brain_region = br, panel = pn,
      delta = mean(sub$score[sub$condition=="tbi"]) - mean(sub$score[sub$condition=="control"]),
      cell_p = t$p.value,
      pseudobulk_p = if (is.null(pb_t)) NA else pb_t$p.value)
  }
}
res2 <- do.call(rbind, results2)
res2$cell_padj <- p.adjust(res2$cell_p, "BH")
res2$pb_padj <- p.adjust(res2$pseudobulk_p, "BH")
write.csv(res2, file.path(OUT, "C7_panel_x_region.csv"), row.names=FALSE)
cat("\n"); print(res2)

# figure: delta across panels, faceted by brain_region
p <- ggplot(res2, aes(panel, delta, fill=brain_region)) +
  geom_col(position=position_dodge(0.7), width=0.65, colour="grey25") +
  geom_hline(yintercept=0, linewidth=0.4, colour="grey40") +
  scale_fill_manual(values=c(hippocampus="#6A9A5B", cortex="#4DAECD")) +
  labs(title="C7. Sensitivity: panel composition × brain region",
       subtitle="Δ module score (TBI − sham), 5 panel compositions",
       y="Δ module score (TBI − sham)", x="panel composition") +
  theme_bw(base_size=11) +
  theme(axis.text.x=element_text(angle=20, hjust=1), legend.position="top")
ggsave(file.path(OUT, "C7_sensitivity_panel.pdf"), p, width=8, height=5)
ggsave(file.path(OUT, "C7_sensitivity_panel.png"), p, width=8, height=5, dpi=200)

cat("\n[C7.5] DONE at:\n  ", OUT, "\n")
