#!/usr/bin/env Rscript
# Devil's Advocate Top 3 - 公共小鼠海马 astrocyte 配体 mRNA 对照
# 目的: 验证人源 GSE209552 的 "ligand mRNA preserved vs communication collapsed" 是否在小鼠成立
# 输入: mouse_integrated_annotated.rds (GSE101901 concussion + GSE180862 FPI)
# 提海马 astrocyte → 看 ligand panel 在 control vs TBI 的 mean expression
suppressMessages({library(Seurat); library(Matrix); library(ggplot2); library(dplyr)})
source("D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/scripts/style_lib_v2.R")

RDS <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/processed_data/integrated/mouse_integrated_annotated.rds"
OUT <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/DA3_public_mouse_ligand"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
LIGANDS <- c("Nrxn1","Nrxn2","Nrxn3","Cadm1","Cadm2","Cadm3","Cntn1","Cntn2","Cntn4","Negr1","Ncam1","Ncam2")

cat("[DA3.1] load + subset hippocampus astrocyte\n")
seu <- readRDS(RDS); DefaultAssay(seu) <- "RNA"
ha <- subset(seu, subset = brain_region == "hippocampus" & celltype_auto == "astrocyte")
rm(seu); gc(verbose=FALSE)
cat("  hippocampal astrocytes:", ncol(ha), "\n")
cat("  by condition:\n"); print(table(ha$condition))

cat("\n[DA3.2] purity gate\n")
expr <- GetAssayData(ha, layer="counts", assay="RNA")
mu <- function(g) if (g %in% rownames(expr)) mean(expr[g,]) else NA
m <- c(Aqp4=mu("Aqp4"), Gfap=mu("Gfap"), Snap25=mu("Snap25"), Syt1=mu("Syt1"))
cat("  marker means:", paste(names(m), round(m,3), sep="=", collapse=", "), "\n")
ratio <- mean(m[c("Aqp4","Gfap")]) / max(0.01, mean(m[c("Snap25","Syt1")]))
cat(sprintf("  astro/neuron ratio = %.2fx (need >=3)\n", ratio))
if (ratio < 3) stop(sprintf("GATE FAIL %.2f", ratio))

cat("\n[DA3.3] per-ligand mean expression by condition\n")
res <- data.frame()
for (g in LIGANDS) {
  if (!g %in% rownames(expr)) next
  for (cond in c("control","tbi")) {
    cells <- which(ha$condition == cond)
    m_c <- mean(expr[g, cells])
    res <- rbind(res, data.frame(gene=g, condition=cond, mean_expr=m_c, n_cells=length(cells)))
  }
}
print(res)

cat("\n[DA3.4] delta + % change\n")
w <- res %>% tidyr::pivot_wider(id_cols=gene, names_from=condition, values_from=mean_expr) %>%
  mutate(delta=tbi-control, pct_change=round(100*delta/abs(control), 1),
         bio_signif = abs(pct_change) >= 15)
print(w)
write.csv(w, file.path(OUT, "DA3_mouse_hippo_astro_ligand.csv"), row.names=FALSE)

cat("\n[DA3.5] sanity: 对比 cell-level Wilcox (per-sample pseudobulk needs animal-level)\n")
ha$condition <- factor(ha$condition, levels=c("control","tbi"))
Idents(ha) <- "condition"
present <- intersect(LIGANDS, rownames(expr))
deg <- FindMarkers(ha, ident.1="tbi", ident.2="control", features=present,
                   logfc.threshold=0, min.pct=0.01, only.pos=FALSE, test.use="wilcox")
deg$gene <- rownames(deg)
deg <- deg[order(deg$gene), ]
print(deg[, c("gene","avg_log2FC","p_val","p_val_adj","pct.1","pct.2")])
write.csv(deg, file.path(OUT, "DA3_DEG_ligands.csv"), row.names=FALSE)

# Figure
cat("\n[DA3.6] figure: ligand mRNA in mouse hippo astrocytes (control vs tbi)\n")
res$condition <- factor(res$condition, levels=c("control","tbi"))
p <- ggplot(res, aes(gene, mean_expr, fill=condition)) +
  geom_col(position=position_dodge(0.78), width=0.7, colour="grey25", linewidth=0.25) +
  geom_text(data=w, aes(x=gene, y=pmax(control, tbi)*1.05,
                        label=paste0(ifelse(pct_change>0,"+",""), pct_change, "%"),
                        fontface=ifelse(bio_signif,"bold","plain")),
            inherit.aes=FALSE, size=2.6) +
  scale_fill_manual(values=c(control=JNI_COL$control, tbi=JNI_COL$tbi)) +
  labs(title="DA3. Mouse hippocampal astrocyte adhesion ligands (mRNA)",
       subtitle="Public mouse integration (GSE101901 + GSE180862, hippocampal astrocyte subset, n>=3 mice/group); % change shown",
       y="mean normalized expression", x=NULL) +
  jni_theme(legend="top")
ggsave(file.path(OUT,"DA3_ligand_bar.pdf"), p, width=8, height=5)
ggsave(file.path(OUT,"DA3_ligand_bar.png"), p, width=8, height=5, dpi=200)

cat("\n[DA3.7] DONE at:\n  ", OUT, "\n")
