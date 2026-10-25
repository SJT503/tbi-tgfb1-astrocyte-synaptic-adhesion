#!/usr/bin/env Rscript
# C5: 神经元亚型黏附塌陷 — 从 v1 mouse_integrated 提 neuron subset, 重聚类找亚型 (兴奋/抑制/海马 CA1/CA3/DG),
#     按亚型分别报黏附模块得分 TBI vs sham, 看哪一亚型塌陷最严重
set.seed(42)
suppressMessages({library(Seurat); library(Matrix); library(ggplot2); library(dplyr)})

RDS <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/processed_data/integrated/mouse_integrated_annotated.rds"
OUT <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/C5_neuron_subtypes"

PANEL <- c("Nrxn1","Nrxn2","Nrxn3","Negr1","Cntn2","Cntn4","Cntn5",
           "Cadm1","Cadm2","Cadm3","Cdh8","Cdh11","Cdh13","Nlgn1","Nlgn2","Nlgn3")

# subtype markers
MARKERS <- list(
  Excitatory_glut     = c("Slc17a7","Slc17a6","Camk2a"),
  Inhibitory_GABA     = c("Gad1","Gad2","Slc32a1"),
  Hippo_CA1_pyramidal = c("Wfs1","Spink8","Pou3f1","Mpped1"),
  Hippo_CA3_pyramidal = c("Prkcd","Cdh24","Calb2","Slc8a2"),
  Hippo_DG_granule    = c("Prox1","C1ql2","Tdo2","Dock10")
)

cat("[C5.1] load + subset neuron\n")
seu <- readRDS(RDS); DefaultAssay(seu) <- "RNA"
neu <- subset(seu, subset = celltype_auto == "neuron")
cat("  neurons:", ncol(neu), " hippocampus:", sum(neu$brain_region == "hippocampus"),
    " cortex:", sum(neu$brain_region == "cortex"), "\n")
rm(seu); gc(verbose=FALSE)

cat("\n[C5.2] AddModuleScore for subtype markers\n")
for (sub in names(MARKERS)) {
  mk <- intersect(MARKERS[[sub]], rownames(neu))
  if (length(mk) < 2) { cat("  skip", sub, ": only", length(mk), "markers\n"); next }
  neu <- AddModuleScore(neu, features=list(mk), name=paste0("SUB_", sub, "_"), seed=42)
}
score_cols <- grep("^SUB_", colnames(neu@meta.data), value=TRUE)
cat("  scored:", paste(score_cols, collapse=","), "\n")

cat("\n[C5.3] assign each cell to highest-scoring subtype\n")
score_mat <- as.matrix(neu@meta.data[, score_cols])
colnames(score_mat) <- gsub("^SUB_|_1$", "", score_cols)
best <- apply(score_mat, 1, function(x) {
  if (max(x) < 0) return("Unassigned")
  colnames(score_mat)[which.max(x)]
})
neu$subtype <- best
cat("  subtype distribution:\n"); print(table(neu$subtype, neu$brain_region))

cat("\n[C5.4] adhesion module score per cell\n")
neu <- AddModuleScore(neu, features=list(intersect(PANEL, rownames(neu))),
                      name="ADH_", seed=42)
neu$ADH_score <- neu$ADH_1

cat("\n[C5.5] per-subtype, per-condition adhesion score (delta TBI - sham)\n")
df <- neu@meta.data %>%
  group_by(subtype, brain_region, condition) %>%
  summarise(n=n(), mean_adh=mean(ADH_score, na.rm=TRUE), sd=sd(ADH_score, na.rm=TRUE),
            .groups="drop")
write.csv(df, file.path(OUT, "C5_per_subtype_adhesion.csv"), row.names=FALSE)
print(df)

# delta
delta <- df %>% tidyr::pivot_wider(id_cols=c(subtype, brain_region),
  names_from=condition, values_from=mean_adh) %>%
  mutate(delta_TBI_vs_sham = tbi - control,
         pct_change = 100 * delta_TBI_vs_sham / abs(control))
write.csv(delta, file.path(OUT, "C5_adhesion_delta.csv"), row.names=FALSE)
cat("\n[C5.6] delta (TBI - sham) per subtype:\n"); print(delta)

# t-test per subtype/region (cell-level)
cat("\n[C5.7] cell-level t-test per subtype × brain_region\n")
test_res <- list()
for (st in unique(neu$subtype)) {
  for (br in unique(neu$brain_region)) {
    cells_st <- which(neu$subtype == st & neu$brain_region == br)
    if (length(cells_st) < 30) next
    sub_meta <- neu@meta.data[cells_st, ]
    if (length(unique(sub_meta$condition)) < 2) next
    t <- tryCatch(t.test(ADH_score ~ condition, data=sub_meta), error=function(e) NULL)
    if (is.null(t)) next
    test_res[[length(test_res)+1]] <- data.frame(
      subtype=st, brain_region=br, n=nrow(sub_meta),
      mean_control=mean(sub_meta$ADH_score[sub_meta$condition=="control"]),
      mean_tbi=mean(sub_meta$ADH_score[sub_meta$condition=="tbi"]),
      t=t$statistic, p_value=t$p.value)
  }
}
test_df <- do.call(rbind, test_res)
test_df$padj_BH <- p.adjust(test_df$p_value, "BH")
test_df <- test_df[order(test_df$p_value), ]
write.csv(test_df, file.path(OUT, "C5_subtype_ttest.csv"), row.names=FALSE)
cat("\nT-test (sorted by p):\n"); print(test_df)

# figure
p1 <- ggplot(df, aes(subtype, mean_adh, fill=condition)) +
  geom_col(position=position_dodge(0.8), width=0.7, colour="grey20") +
  facet_wrap(~brain_region) +
  scale_fill_manual(values=c(control="#4DAECD", tbi="#E07B54")) +
  labs(title="C5. Neuronal adhesion module score by subtype",
       subtitle="Public hippocampal+cortical mouse TBI (mouse_integrated_annotated.rds)",
       y="Adhesion module score (16 genes)", x=NULL) +
  theme_bw(base_size=10) +
  theme(axis.text.x=element_text(angle=30, hjust=1), legend.position="top")
ggsave(file.path(OUT, "C5_subtype_adhesion.pdf"), p1, width=9, height=5)
ggsave(file.path(OUT, "C5_subtype_adhesion.png"), p1, width=9, height=5, dpi=200)

cat("\n[C5.8] DONE at:\n  ", OUT, "\n")
