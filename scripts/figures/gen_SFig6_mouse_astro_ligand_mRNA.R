#!/usr/bin/env Rscript
# SFig 6 — Mouse hippocampal astrocyte adhesion ligands (mRNA) from saved DA3 tables
suppressMessages({library(ggplot2); library(dplyr); library(tidyr)})
source("D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/scripts/style_lib_v2.R")

BASE <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI"
IN   <- file.path(BASE, "analysis/DA3_public_mouse_ligand")
OUT  <- file.path(BASE, "supplementary")

w   <- read.csv(file.path(IN, "DA3_mouse_hippo_astro_ligand.csv"))
deg <- read.csv(file.path(IN, "DA3_DEG_ligands.csv"))

w <- merge(w, deg[, c("gene", "p_val_adj")], by="gene", all.x=TRUE)
w$sig <- ifelse(is.na(w$p_val_adj), "",
          ifelse(w$p_val_adj < 0.001, "***",
          ifelse(w$p_val_adj < 0.01,  "**",
          ifelse(w$p_val_adj < 0.05,  "*", ""))))

gene_order <- c("Nrxn1","Nrxn2","Nrxn3","Cadm1","Cadm2","Cadm3","Cntn1","Cntn2","Cntn4","Negr1","Ncam1","Ncam2")
w$gene <- factor(w$gene, levels=gene_order)

res <- w %>%
  select(gene, control, tbi) %>%
  pivot_longer(cols=c("control","tbi"), names_to="condition", values_to="mean_expr") %>%
  mutate(condition_lab = factor(ifelse(condition=="control","Sham","TBI"), levels=c("Sham","TBI")))

max_y <- max(w$control, w$tbi, na.rm=TRUE)
w$label <- paste0(ifelse(w$pct_change>0, "+", ""), w$pct_change, "%",
                  ifelse(w$sig=="", "", paste0("\n", w$sig)))

p <- ggplot(res, aes(gene, mean_expr, fill=condition_lab)) +
  geom_col(position=position_dodge(0.76), width=0.68, colour="grey25", linewidth=0.25) +
  geom_text(data=w, aes(x=gene, y=pmax(control, tbi) + max_y * 0.055,
                        label=label,
                        fontface=ifelse(bio_signif, "bold", "plain")),
            inherit.aes=FALSE, size=2.6, lineheight=0.85) +
  scale_fill_manual(values=c(Sham=JNI_COL$control, TBI=JNI_COL$tbi)) +
  scale_y_continuous(limits=c(0, max_y * 1.18), expand=expansion(mult=c(0, 0.02))) +
  labs(title="Supplementary Figure S6. Astrocytic adhesion-ligand mRNA in mouse hippocampus after TBI",
       subtitle="Public mouse hippocampal astrocytes (GSE101901 + GSE180862, n = 9 Sham vs 13 TBI); bold = |delta%| >= 15%; *** padj < 0.001",
       y="Mean normalized expression", x=NULL) +
  jni_theme(legend="top") +
  theme(plot.title=element_text(face="bold", size=10.5),
        plot.subtitle=element_text(size=8, colour="grey35"),
        plot.margin=margin(4, 6, 4, 4, "mm"))

ggsave(file.path(OUT, "FigS6_mouse_astro_ligand_mRNA.pdf"), p, width=195, height=120, units="mm")
ggsave(file.path(OUT, "FigS6_mouse_astro_ligand_mRNA.png"), p, width=195, height=120, units="mm", dpi=300)
cat("Saved FigS6_mouse_astro_ligand_mRNA.pdf/.png\n")
