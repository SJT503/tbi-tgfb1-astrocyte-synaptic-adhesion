#!/usr/bin/env Rscript
# Step 5 headline figure: TGFBR1 transcript by cell type (mouse + human)
suppressMessages({library(ggplot2)})
OUT <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/step5_out"

# --- mouse (from 5a output) ---
m <- read.csv(file.path(OUT,"5a_mouse_tgfb_byCelltype.csv"))
m <- m[m$gene=="Tgfbr1",]
m$celltype <- factor(c(astrocyte="Astrocyte",neuron="Neuron",
  microglia_macrophage="Microglia",oligodendrocyte="Oligo",
  endothelial="Endo",opc="OPC")[m$celltype],
  levels=c("Microglia","OPC","Endo","Oligo","Astrocyte","Neuron"))
m$species <- "Mouse (GSE101901+180862)"

# --- human (precomputed CSV) ---
hpath <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/wet_lab_data/results/figures/Fig_Human_TGFb/human_tgfb_byCelltype.csv"
h <- read.csv(hpath)
h <- h[h$gene=="TGFBR1",]
h$group <- ifelse(h$group=="CTL","control","tbi")
h$celltype <- factor(c(Astro="Astrocyte",Neuron="Neuron",Micro="Microglia",
  Oligo="Oligo",Endo="Endo")[h$celltype],
  levels=c("Microglia","OPC","Endo","Oligo","Astrocyte","Neuron"))
h$species <- "Human (GSE209552)"
h$det_pct <- h$det

d <- rbind(m[,c("celltype","group","mean","det_pct","species")],
           h[,c("celltype","group","mean","det_pct","species")])
d <- d[!is.na(d$celltype),]

p <- ggplot(d, aes(celltype, mean, fill=group)) +
  geom_col(position=position_dodge(0.75), width=0.68, colour="grey20", linewidth=0.25) +
  facet_wrap(~species, scales="free_x") +
  scale_fill_manual(values=c(control="#4DAECD", tbi="#E07B54"), name=NULL) +
  coord_flip() +
  labs(title="TGFBR1 transcript is microglia-predominant in both species",
       subtitle="single-cell / single-nucleus RNA; mean normalized expression",
       x=NULL, y="Mean expression") +
  theme_bw(base_size=9) +
  theme(plot.title=element_text(face="bold", size=9.5),
        plot.subtitle=element_text(size=7.6, colour="grey35"),
        strip.background=element_rect(fill="grey92"),
        strip.text=element_text(face="bold", size=7.8),
        legend.position="top", legend.key.size=unit(3.2,"mm"),
        panel.grid.minor=element_blank())

ggsave(file.path(OUT,"step5_TGFBR1_byCelltype.pdf"), p, width=150, height=70, units="mm")
ggsave(file.path(OUT,"step5_TGFBR1_byCelltype.png"), p, width=150, height=70, units="mm", dpi=300)
cat("saved step5_TGFBR1_byCelltype\n")
print(d[order(d$species,d$celltype),], row.names=FALSE)
