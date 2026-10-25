#!/usr/bin/env Rscript
# Candidate 1 — within astrocytes, is TGF-β/Smad signalling activity
# anti-correlated with synaptic-adhesion-molecule expression?
# Tests the 6a upstream link: astrocytic TGF-βRI ON -> adhesion ligand DOWN.
# Real processed data only.
suppressMessages({library(Seurat); library(ggplot2)})
PD  <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/processed_data"
OUT <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/cand_out"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

TGFB <- c("Serpine1","Smad7","Pmepa1","Skil","Junb","Tgfbi","Bambi","Id1","Id2","Id3","Ltbp3","Smad3")
ADH  <- c("Ncam1","Nrxn1","Nrxn3","Cntn1","Negr1","Cadm1")

mo <- readRDS(file.path(PD,"integrated/mouse_integrated_annotated.rds"))
DefaultAssay(mo) <- "RNA"
ast <- subset(mo, subset = celltype_auto == "astrocyte")
m <- GetAssayData(ast, assay="RNA", layer="data")
r <- mean(m["Aqp4",]) / (mean(m["Snap25",])+1e-6)
cat(sprintf("PURITY astrocyte Aqp4/Snap25 = %.1f %s\n", r, ifelse(r>=3,"PASS","** FAIL **")))
if(r < 3) stop("PURITY FAIL: astrocyte subset impure")
cat("astrocytes:", ncol(ast), " control/tbi:", paste(table(ast$condition),collapse="/"), "\n")

tg <- TGFB[TGFB %in% rownames(ast)]; ad <- ADH[ADH %in% rownames(ast)]
cat("TGF-β signalling genes used:", paste(tg,collapse=","), "\n")
cat("adhesion genes used:", paste(ad,collapse=","), "\n")
set.seed(1)
ast <- AddModuleScore(ast, features=list(TGFb=tg, Adh=ad), name="M_", assay="RNA")
df <- data.frame(tgfb=ast$M_1, adh=ast$M_2, condition=ast$condition, dataset=ast$dataset)

# correlations: overall + within each condition (controls for TBI shifting both)
res <- list()
for(grp in c("ALL","control","tbi")){
  sub <- if(grp=="ALL") df else df[df$condition==grp,]
  ct <- cor.test(sub$tgfb, sub$adh, method="spearman")
  res[[grp]] <- data.frame(group=grp, n=nrow(sub),
    spearman_rho=round(unname(ct$estimate),4), p=signif(ct$p.value,3))
}
res <- do.call(rbind, res)
write.csv(res, file.path(OUT,"cand1_correlations.csv"), row.names=FALSE)
write.csv(df,  file.path(OUT,"cand1_astrocyte_scores.csv"), row.names=FALSE)
cat("\n=== TGF-β signalling vs adhesion module, astrocyte-level Spearman ===\n")
print(res, row.names=FALSE)

# figure
p <- ggplot(df, aes(tgfb, adh)) +
  geom_point(aes(colour=condition), size=.35, alpha=.3) +
  geom_smooth(method="lm", se=TRUE, colour="#b8412e", linewidth=.7) +
  scale_colour_manual(values=c(control="#4DAECD", tbi="#E07B54")) +
  labs(title="Astrocyte TGF-β signalling vs synaptic-adhesion module",
       subtitle=sprintf("mouse public atlas, %d astrocytes; Spearman rho(all)=%.3f, p=%.2g",
                        nrow(df), res$spearman_rho[res$group=="ALL"], res$p[res$group=="ALL"]),
       x="TGF-β / Smad signalling module score", y="Synaptic-adhesion module score") +
  theme_bw(base_size=9) +
  theme(plot.title=element_text(face="bold",size=9.5),
        plot.subtitle=element_text(size=7.4,colour="grey35"),
        legend.position="top", legend.title=element_blank(),
        panel.grid.minor=element_blank())
ggsave(file.path(OUT,"cand1_tgfb_vs_adhesion.pdf"), p, width=110, height=105, units="mm")
ggsave(file.path(OUT,"cand1_tgfb_vs_adhesion.png"), p, width=110, height=105, units="mm", dpi=300)
cat("\nsaved cand1 figure\nDONE cand1\n")
