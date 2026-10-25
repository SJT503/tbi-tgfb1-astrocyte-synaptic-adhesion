#!/usr/bin/env Rscript
# Step 6 supporting analysis — neuronal survival vs intrinsic-apoptosis vs adhesion
# module scores, control vs TBI (public mouse atlas). Transcriptomic-layer support
# for the anoikis arm. Real processed data only.
suppressMessages({library(Seurat); library(ggplot2)})
PD  <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/processed_data"
OUT <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/step6_out"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

mo <- readRDS(file.path(PD,"integrated/mouse_integrated_annotated.rds"))
DefaultAssay(mo) <- "RNA"

# --- subset neurons + purity gate ---
neu <- subset(mo, subset = celltype_auto == "neuron")
m <- GetAssayData(neu, assay="RNA", layer="data")
r <- mean(m["Snap25",]) / (mean(m["Aqp4",])+1e-6)
cat(sprintf("PURITY neuron Snap25/Aqp4 = %.1f %s\n", r, ifelse(r>=3,"PASS","** FAIL **")))
if(r < 3) stop("PURITY FAIL: neuron subset not pure")
cat("neurons:", ncol(neu), " control/tbi:", paste(table(neu$condition),collapse="/"), "\n")

# --- gene modules (mouse symbols) ---
mods <- list(
  IntrinsicApoptosis = c("Bax","Bak1","Bbc3","Pmaip1","Apaf1","Casp9","Casp3","Cycs"),
  Antiapoptotic      = c("Bcl2","Bcl2l1","Mcl1","Bcl2l2","Xiap"),
  Adhesion           = c("Ncam1","Nrxn1","Nrxn3","Cntn1","Negr1","Cadm1"))
for(nm in names(mods)){
  g <- mods[[nm]]; present <- g[g %in% rownames(neu)]
  cat(sprintf("  %s: %d/%d genes -> %s\n", nm, length(present), length(g),
              paste(present,collapse=",")))
  mods[[nm]] <- present
}

set.seed(1)
neu <- AddModuleScore(neu, features=mods, name="MOD_", assay="RNA")
sc <- grep("^MOD_", colnames(neu@meta.data), value=TRUE)
names(sc) <- names(mods)

# --- per-cell test, by dataset ---
res <- list()
for(ds in unique(neu$dataset)){
  for(nm in names(mods)){
    col <- sc[nm]
    sub <- neu@meta.data[neu$dataset==ds,]
    ctl <- sub[[col]][sub$condition=="control"]; tbi <- sub[[col]][sub$condition=="tbi"]
    if(length(ctl)<10 || length(tbi)<10) next
    w <- wilcox.test(ctl, tbi)
    res[[length(res)+1]] <- data.frame(dataset=ds, module=nm,
      n_ctl=length(ctl), n_tbi=length(tbi),
      mean_ctl=round(mean(ctl),4), mean_tbi=round(mean(tbi),4),
      delta=round(mean(tbi)-mean(ctl),4), wilcox_p=signif(w$p.value,3))
  }
}
res <- do.call(rbind, res)
write.csv(res, file.path(OUT,"step6_neuron_module_percell.csv"), row.names=FALSE)
cat("\n=== per-cell module score: control vs TBI ===\n"); print(res, row.names=FALSE)

# --- pseudobulk per sample (GSE101901 has Sham1-3/TBI1-3) ---
pb <- aggregate(neu@meta.data[,sc], by=list(sample=neu$orig.ident,
                cond=neu$condition, dataset=neu$dataset), FUN=mean)
colnames(pb)[match(sc,colnames(pb))] <- names(mods)
write.csv(pb, file.path(OUT,"step6_neuron_module_pseudobulk.csv"), row.names=FALSE)
cat("\n=== pseudobulk (mean module score per sample) ===\n"); print(pb, row.names=FALSE)
g101 <- pb[pb$dataset=="GSE101901",]
if(nrow(g101)>=4){
  cat("\nGSE101901 pseudobulk t-test (Sham n=3 vs TBI n=3):\n")
  for(nm in names(mods)){
    c0 <- g101[[nm]][g101$cond=="control"]; c1 <- g101[[nm]][g101$cond=="tbi"]
    if(length(c0)>=2 && length(c1)>=2){
      tt <- t.test(c0,c1)
      cat(sprintf("  %s: ctl=%.4f tbi=%.4f  p=%.3f\n", nm, mean(c0), mean(c1), tt$p.value))
    }
  }
}

# --- figure ---
df <- do.call(rbind, lapply(names(mods), function(nm)
  data.frame(score=neu@meta.data[[sc[nm]]], module=nm,
             condition=neu$condition, dataset=neu$dataset)))
df$module <- factor(df$module, levels=c("Adhesion","Antiapoptotic","IntrinsicApoptosis"))
p <- ggplot(df, aes(condition, score, fill=condition)) +
  geom_violin(scale="width", trim=TRUE, linewidth=0.3) +
  geom_boxplot(width=0.16, outlier.shape=NA, fill="white", linewidth=0.3) +
  facet_grid(dataset~module, scales="free_y") +
  scale_fill_manual(values=c(control="#4DAECD", tbi="#E07B54")) +
  labs(title="Neuronal module scores: control vs TBI (public mouse atlas)",
       subtitle="transcriptomic-layer support for the anoikis arm",
       x=NULL, y="Module score") +
  theme_bw(base_size=9) +
  theme(plot.title=element_text(face="bold",size=9.5),
        plot.subtitle=element_text(size=7.6,colour="grey35"),
        legend.position="none", strip.background=element_rect(fill="grey92"),
        strip.text=element_text(face="bold",size=7.6), panel.grid.minor=element_blank())
ggsave(file.path(OUT,"step6_neuron_modules.pdf"), p, width=150, height=120, units="mm")
ggsave(file.path(OUT,"step6_neuron_modules.png"), p, width=150, height=120, units="mm", dpi=300)
cat("\nsaved step6_neuron_modules figure\nDONE\n")
