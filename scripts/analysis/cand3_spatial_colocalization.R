#!/usr/bin/env Rscript
# Candidate 3 — spatial co-localization: do tissue regions with high astrocyte
# TGF-β signalling show LOW neuronal synaptic-adhesion expression?
# GSE226208 Visium. Real processed data only. Spot-resolution caveat noted.
suppressMessages({library(Seurat); library(ggplot2)})
PD  <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/processed_data"
OUT <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/cand_out"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

TGFB   <- c("Serpine1","Smad7","Pmepa1","Skil","Junb","Tgfbi","Bambi","Id1","Id2","Id3","Ltbp3","Smad3")
ADH    <- c("Ncam1","Nrxn1","Nrxn3","Cntn1","Negr1","Cadm1")
ASTRO  <- c("Aqp4","Gfap","Slc1a3","Aldh1l1")
NEURON <- c("Snap25","Rbfox3","Syt1")

partial_cor <- function(x, y, z){            # partial Spearman of x~y given z
  rx <- residuals(lm(rank(x)~rank(z))); ry <- residuals(lm(rank(y)~rank(z)))
  ct <- cor.test(rx, ry, method="spearman"); c(rho=unname(ct$estimate), p=ct$p.value)
}

analyse <- function(path, tag){
  if(!file.exists(path)){ cat(tag,"missing\n"); return(NULL) }
  sp <- readRDS(path); DefaultAssay(sp) <- "Spatial"
  sp <- NormalizeData(sp, verbose=FALSE)
  m <- GetAssayData(sp, assay="Spatial", layer="data")
  sc <- function(gs){ g<-gs[gs%in%rownames(m)]; if(length(g)) Matrix::colMeans(m[g,,drop=FALSE]) else rep(NA,ncol(m)) }
  d <- data.frame(tgfb=sc(TGFB), adh=sc(ADH), astro=sc(ASTRO), neuron=sc(NEURON))
  cat(sprintf("\n[%s] %d spots\n", tag, nrow(d)))
  cm <- round(cor(d, method="spearman", use="pairwise.complete.obs"), 3)
  cat("Spearman matrix:\n"); print(cm)
  # key test: TGF-β signalling vs neuronal adhesion, controlling for neuron content
  pc <- partial_cor(d$tgfb, d$adh, d$neuron)
  cat(sprintf("PARTIAL Spearman  tgfb ~ adh | neuron : rho=%.3f  p=%.2g\n", pc["rho"], pc["p"]))
  pc2 <- partial_cor(d$astro, d$adh, d$neuron)
  cat(sprintf("PARTIAL Spearman  astro ~ adh | neuron: rho=%.3f  p=%.2g\n", pc2["rho"], pc2["p"]))
  d$region <- tag
  list(scores=d, cormat=cm,
       summary=data.frame(region=tag, n=nrow(d),
         tgfb_adh_raw=cm["tgfb","adh"],
         tgfb_adh_partial=round(pc["rho"],3), tgfb_adh_partial_p=signif(pc["p"],3),
         astro_adh_partial=round(pc2["rho"],3)))
}

r1 <- analyse(file.path(PD,"spatial/visium_tbi_3dpi.rds"), "TBI_3dpi")
r0 <- analyse(file.path(PD,"spatial/visium_intact.rds"),   "intact")

alld  <- do.call(rbind, lapply(list(r0,r1), function(x) if(!is.null(x)) x$scores))
summ  <- do.call(rbind, lapply(list(r0,r1), function(x) if(!is.null(x)) x$summary))
write.csv(alld, file.path(OUT,"cand3_spatial_spot_scores.csv"), row.names=FALSE)
write.csv(summ, file.path(OUT,"cand3_spatial_summary.csv"), row.names=FALSE)
cat("\n=== summary ===\n"); print(summ, row.names=FALSE)

if(!is.null(alld)){
  p <- ggplot(alld, aes(tgfb, adh)) +
    geom_point(size=.5, alpha=.35, colour="#E07B54") +
    geom_smooth(method="lm", se=TRUE, colour="#1c3553", linewidth=.7) +
    facet_wrap(~region, scales="free") +
    labs(title="Spatial: astrocyte TGF-β signalling vs neuronal adhesion (per spot)",
         subtitle="GSE226208 Visium; spot resolution = multi-cell — see partial correlations in CSV",
         x="TGF-β / Smad signalling (spot)", y="Synaptic-adhesion module (spot)") +
    theme_bw(base_size=9) +
    theme(plot.title=element_text(face="bold",size=9.3),
          plot.subtitle=element_text(size=7.2,colour="grey35"),
          strip.background=element_rect(fill="grey92"), panel.grid.minor=element_blank())
  ggsave(file.path(OUT,"cand3_spatial_tgfb_vs_adhesion.pdf"), p, width=150, height=80, units="mm")
  ggsave(file.path(OUT,"cand3_spatial_tgfb_vs_adhesion.png"), p, width=150, height=80, units="mm", dpi=300)
  cat("saved cand3 figure\n")
}
cat("DONE cand3\n")
