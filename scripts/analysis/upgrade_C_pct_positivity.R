#!/usr/bin/env Rscript
# upgrade_C_pct_positivity.R — v2.2 升级 C
# % positivity report on Astrocyte & Neuron, three groups (descriptive, no p-value)
# 输入: all_cells_annotated.rds
# 输出: % positive cells per gene per celltype per group + 柱状图
suppressMessages({library(Seurat); library(Matrix); library(ggplot2); library(dplyr)})

RDS <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/wet_lab_data/results/seurat/all_cells_annotated.rds"
OUT <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/upgrade_C_pct_positivity"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

PANEL_NEURON <- c("Nrxn1","Nrxn2","Nrxn3","Negr1","Cntn2","Cntn4","Cntn5",
                  "Cadm1","Cadm2","Cadm3","Cdh8","Cdh11","Cdh13",
                  "Nlgn1","Nlgn2","Nlgn3")
PANEL_ASTRO <- c("Nrxn1","Cadm1","Cadm2","Nlgn1","Nlgn3","Cntn1","Cntn2")  # astrocytic adhesion partners

cat("[C.1] load\n")
seu <- readRDS(RDS); DefaultAssay(seu) <- "RNA"
expr <- GetAssayData(seu, layer="counts", assay="RNA")
grp <- as.character(seu$stim); ct <- as.character(seu$celltype)
grp_levels <- c("Sham","CCI_Scr","CCI_TGFbetaR1_KO")

pct_pos <- function(genes, cells_idx) {
  out <- sapply(genes, function(g) {
    if (!g %in% rownames(expr)) return(NA_real_)
    100 * mean(expr[g, cells_idx] >= 1)
  })
  out
}

res <- list()
for (ctype in c("Neuron","Astrocyte")) {
  panel <- if (ctype == "Neuron") PANEL_NEURON else PANEL_ASTRO
  for (g in grp_levels) {
    cells_idx <- which(ct == ctype & grp == g)
    if (length(cells_idx) < 5) next
    v <- pct_pos(panel, cells_idx)
    res[[length(res)+1]] <- data.frame(
      celltype = ctype, group = g, gene = names(v),
      pct = unname(v), n_cells = length(cells_idx))
  }
}
df <- do.call(rbind, res)
df$group <- factor(df$group, levels=grp_levels)
write.csv(df, file.path(OUT, "C_pct_positivity.csv"), row.names=FALSE)
print(head(df, 20))

# 关键三组比较: 每基因 KO vs CCI_Scr 恢复方向
wide <- df %>% tidyr::pivot_wider(names_from=group, values_from=pct,
                                  id_cols=c(celltype, gene))
wide$delta_TBI_vs_Sham <- wide$CCI_Scr - wide$Sham
wide$delta_KO_recovery <- wide$CCI_TGFbetaR1_KO - wide$CCI_Scr
wide$pct_rescue <- wide$delta_KO_recovery / (wide$Sham - wide$CCI_Scr) * 100
write.csv(wide, file.path(OUT, "C_recovery_summary.csv"), row.names=FALSE)
cat("\n[C.2] recovery summary\n")
print(wide)

# 图: per-celltype facet, 三组柱
cat("\n[C.3] bar plot\n")
p <- ggplot(df, aes(gene, pct, fill=group)) +
  geom_col(position=position_dodge(0.8), width=0.7, colour="grey25", linewidth=0.2) +
  facet_wrap(~celltype, scales="free", ncol=1) +
  scale_fill_manual(values=c(Sham="#4DAECD", CCI_Scr="#E07B54", CCI_TGFbetaR1_KO="#6A9A5B")) +
  labs(title="C. % positive cells across three conditions (descriptive)",
       subtitle="≥1 UMI = positive; n=1 pooled library per group, no p-value reported",
       x=NULL, y="% positive cells") +
  theme_bw(base_size=10) +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        legend.position="top")
ggsave(file.path(OUT, "C_pct_positivity_bars.pdf"), p, width=8, height=6)
ggsave(file.path(OUT, "C_pct_positivity_bars.png"), p, width=8, height=6, dpi=200)

cat("\n[C.4] DONE. Outputs at:\n  ", OUT, "\n")
