#!/usr/bin/env Rscript
# DA3 镜像版 - Mouse CORTICAL astrocyte ligand mRNA (对照海马 DA3)
# 目的: 在小鼠内部直接对比 hippocampus vs cortex 区域差异
# 输入: mouse_integrated_annotated.rds, subset cortex + astrocyte
set.seed(42)
suppressMessages({library(Seurat); library(Matrix); library(ggplot2); library(dplyr)})
source("D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/scripts/style_lib_v2.R")

RDS <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/processed_data/integrated/mouse_integrated_annotated.rds"
OUT <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/DA3c_mouse_cortex_ligand"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
LIGANDS <- c("Nrxn1","Nrxn2","Nrxn3","Cadm1","Cadm2","Cadm3","Cntn1","Cntn2","Cntn4","Negr1","Ncam1","Ncam2")

cat("[C.1] load + subset CORTEX astrocyte\n")
seu <- readRDS(RDS); DefaultAssay(seu) <- "RNA"
ca <- subset(seu, subset = brain_region == "cortex" & celltype_auto == "astrocyte")
rm(seu); gc(verbose=FALSE)
cat("  cortical astrocytes:", ncol(ca), "\n")
cat("  by condition:\n"); print(table(ca$condition))
cat("  by animal x condition:\n"); print(table(ca$animal, ca$condition))

cat("\n[C.2] purity gate\n")
expr <- GetAssayData(ca, layer="counts", assay="RNA")
mu <- function(g) if (g %in% rownames(expr)) mean(expr[g,]) else NA
m <- c(Aqp4=mu("Aqp4"), Gfap=mu("Gfap"), Snap25=mu("Snap25"), Syt1=mu("Syt1"))
cat("  marker means:", paste(names(m), round(m,3), sep="=", collapse=", "), "\n")
ratio <- mean(m[c("Aqp4","Gfap")]) / max(0.01, mean(m[c("Snap25","Syt1")]))
cat(sprintf("  astro/neuron ratio = %.2fx (need >=3)\n", ratio))
if (ratio < 3) stop(sprintf("GATE FAIL %.2f", ratio))

cat("\n[C.3] per-ligand mean expression by condition\n")
res <- data.frame()
for (g in LIGANDS) {
  if (!g %in% rownames(expr)) next
  for (cond in c("control","tbi")) {
    cells <- which(ca$condition == cond)
    res <- rbind(res, data.frame(gene=g, condition=cond, mean_expr=mean(expr[g, cells]), n_cells=length(cells)))
  }
}
print(res)

cat("\n[C.4] delta + % change\n")
w <- res %>% tidyr::pivot_wider(id_cols=gene, names_from=condition, values_from=mean_expr) %>%
  mutate(delta=tbi-control, pct_change=round(100*delta/abs(control), 1),
         bio_signif = abs(pct_change) >= 15)
print(w)
write.csv(w, file.path(OUT, "DA3c_mouse_cortex_astro_ligand.csv"), row.names=FALSE)

cat("\n[C.5] Wilcox per-ligand\n")
ca$condition <- factor(ca$condition, levels=c("control","tbi"))
Idents(ca) <- "condition"
present <- intersect(LIGANDS, rownames(expr))
deg <- FindMarkers(ca, ident.1="tbi", ident.2="control", features=present,
                   logfc.threshold=0, min.pct=0.01, only.pos=FALSE, test.use="wilcox")
deg$gene <- rownames(deg); deg <- deg[order(deg$gene), ]
print(deg[, c("gene","avg_log2FC","p_val","p_val_adj","pct.1","pct.2")])
write.csv(deg, file.path(OUT, "DA3c_DEG_ligands.csv"), row.names=FALSE)

cat("\n[C.6] figure\n")
res$condition <- factor(res$condition, levels=c("control","tbi"))
p <- ggplot(res, aes(gene, mean_expr, fill=condition)) +
  geom_col(position=position_dodge(0.78), width=0.7, colour="grey25", linewidth=0.25) +
  geom_text(data=w, aes(x=gene, y=pmax(control, tbi)*1.05,
                        label=paste0(ifelse(pct_change>0,"+",""), pct_change, "%"),
                        fontface=ifelse(bio_signif,"bold","plain")),
            inherit.aes=FALSE, size=2.6) +
  scale_fill_manual(values=c(control=JNI_COL$control, tbi=JNI_COL$tbi)) +
  labs(title="DA3c. Mouse CORTICAL astrocyte adhesion ligands (mRNA)",
       subtitle="Mirror control of DA3 hippocampus; GSE180862 cortex subset; same ligand panel",
       y="mean normalized expression", x=NULL) +
  jni_theme(legend="top")
ggsave(file.path(OUT,"DA3c_ligand_bar.pdf"), p, width=8, height=5)
ggsave(file.path(OUT,"DA3c_ligand_bar.png"), p, width=8, height=5, dpi=200)

cat("\n[C.7] DONE at:\n  ", OUT, "\n")
