#!/usr/bin/env Rscript
# Supplementary figure — computational support for the astrocyte-neuron adhesion axis
#  A  TGFBR1 transcript by cell type (mouse + human)
#  B  astrocyte->neuron synaptic-adhesion communication, control vs TBI (public CellChat)
#  C  astrocyte->DG-neuron adhesion-pathway strength, Sham/CCI/KO (self snRNA, descriptive)
#  D  neuronal adhesion module score change across TBI severity gradient
suppressMessages({library(ggplot2); library(patchwork)})
BASE <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project"
S5  <- file.path(BASE,"JNI/analysis/step5_out")
S6  <- file.path(BASE,"JNI/analysis/step6_out")
OUT <- file.path(BASE,"JNI/supplementary")
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
col_ctl <- "#4DAECD"; col_tbi <- "#E07B54"; col_ko <- "#6A9A5B"

pub <- theme_bw(base_size=9) +
  theme(plot.title=element_text(face="bold",size=9.5),
        plot.subtitle=element_text(size=7.2,colour="grey35"),
        strip.background=element_rect(fill="grey92"),
        strip.text=element_text(face="bold",size=7.6),
        panel.grid.minor=element_blank(),
        legend.key.size=unit(3.2,"mm"), legend.title=element_blank())

# ---- Panel A : TGFBR1 by cell type ----
m <- read.csv(file.path(S5,"5a_mouse_tgfb_byCelltype.csv")); m <- m[m$gene=="Tgfbr1",]
m$celltype <- c(astrocyte="Astrocyte",neuron="Neuron",microglia_macrophage="Microglia",
  oligodendrocyte="Oligo",endothelial="Endo",opc="OPC")[m$celltype]
m$species <- "Mouse"
h <- read.csv(file.path(BASE,"wet_lab_data/results/figures/Fig_Human_TGFb/human_tgfb_byCelltype.csv"))
h <- h[h$gene=="TGFBR1",]; h$group <- ifelse(h$group=="CTL","control","tbi")
h$celltype <- c(Astro="Astrocyte",Neuron="Neuron",Micro="Microglia",Oligo="Oligo",Endo="Endo")[h$celltype]
h$species <- "Human"
A <- rbind(m[,c("celltype","group","mean","species")], h[,c("celltype","group","mean","species")])
A$celltype <- factor(A$celltype, levels=c("Microglia","OPC","Endo","Oligo","Astrocyte","Neuron"))
A$species  <- factor(A$species, levels=c("Mouse","Human"))
A <- A[!is.na(A$celltype),]
pA <- ggplot(A, aes(celltype, mean, fill=group)) +
  geom_col(position=position_dodge(0.74), width=0.66, colour="grey25", linewidth=0.25) +
  facet_wrap(~species, scales="free_x") + coord_flip() +
  scale_fill_manual(values=c(control=col_ctl, tbi=col_tbi)) +
  labs(title="A  TGFBR1 transcript by cell type",
       subtitle="single-cell/nucleus RNA; receptor mRNA is microglia-predominant in both species",
       x=NULL, y="Mean expression") + pub + theme(legend.position="top")

# ---- Panel B : public CellChat astrocyte->neuron adhesion pathways ----
b <- read.csv(file.path(S5,"5b_astro_to_neuron_pathways.csv"))
adh <- c("NRXN","NCAM","NEGR","NGL","CNTN","CADM")
b <- b[b$pathway_name %in% adh,]
b$pathway_name <- factor(b$pathway_name, levels=adh)
b$dir <- ifelse(b$log2FC<0,"down","up")
pB <- ggplot(b, aes(pathway_name, log2FC, fill=dir)) +
  geom_col(width=0.62, colour="grey25", linewidth=0.25) +
  geom_hline(yintercept=0, linewidth=0.4, colour="grey40") +
  scale_fill_manual(values=c(down=col_ctl, up=col_tbi), guide="none") +
  labs(title="B  Astrocyte->neuron adhesion: mouse CellChat",
       subtitle="public mouse, TBI vs control; NRXN/NEGR/CNTN/NGL collapse (NCAM not)",
       x=NULL, y="log2(TBI / control)") + pub

# ---- Panel E : human CellChat astrocyte->neuron adhesion (strict purity) ----
e <- read.csv(file.path(BASE,"JNI/analysis/cand_out/cand2b_human_astro2neuron_pathways.csv"))
e <- e[e$adhesion==TRUE & !is.na(e$pct_change) & is.finite(e$pct_change),]
e$pathway <- factor(e$pathway, levels=e$pathway[order(e$pct_change)])
e$dir <- ifelse(e$pct_change<0,"down","up")
pE <- ggplot(e, aes(pathway, pct_change, fill=dir)) +
  geom_col(width=0.62, colour="grey25", linewidth=0.25) +
  geom_hline(yintercept=0, linewidth=0.4, colour="grey40") +
  scale_fill_manual(values=c(down=col_ctl, up=col_tbi), guide="none") +
  labs(title="E  Astrocyte->neuron adhesion: human CellChat",
       subtitle="human GSE209552, TBI vs control; NRXN-led collapse, NCAM discordant",
       x=NULL, y="% change (TBI vs control)") + pub

# ---- Panel C : self snRNA 3-group adhesion pathway strength ----
cc <- read.csv(file.path(BASE,"wet_lab_data/results/cellchat/Ast_to_DG_pathway_strength_3groups.csv"))
cc <- cc[,c("Pathway","Sham","CCI_Scr","KO")]
cl <- reshape(cc, varying=c("Sham","CCI_Scr","KO"), v.names="strength",
              timevar="group", times=c("Sham","CCI_Scr","KO"), direction="long")
cl$group <- factor(cl$group, levels=c("Sham","CCI_Scr","KO"))
cl$Pathway <- factor(cl$Pathway, levels=cc$Pathway[order(-cc$Sham)])
pC <- ggplot(cl, aes(Pathway, strength, fill=group)) +
  geom_col(position=position_dodge(0.76), width=0.7, colour="grey25", linewidth=0.25) +
  scale_fill_manual(values=c(Sham=col_ctl, CCI_Scr=col_tbi, KO=col_ko)) +
  labs(title="C  Astrocyte->neuron adhesion: self snRNA (3 groups)",
       subtitle="descriptive only (n=1 pooled library/group); CCI lowers, astrocytic TGFBRI-KO restores",
       x=NULL, y="Communication strength") + pub + theme(legend.position="top")

# ---- Panel D : neuronal adhesion module across severity gradient ----
d <- read.csv(file.path(S6,"step6b_severity_gradient.csv"))
d <- d[d$module=="Adhesion",]
lab <- c(GSE290150="Moderate",GSE271769="Severe-1",GSE277487="Severe-2",GSE209552="Human")
d$dataset <- factor(lab[d$dataset], levels=lab)
d$dir <- ifelse(d$delta<0,"down","up")
d$sig <- ifelse(d$pseudobulk_p<0.001,"***",ifelse(d$pseudobulk_p<0.01,"**",
          ifelse(d$pseudobulk_p<0.05,"*","ns")))
pD <- ggplot(d, aes(dataset, delta, fill=dir)) +
  geom_col(width=0.6, colour="grey25", linewidth=0.25) +
  geom_hline(yintercept=0, linewidth=0.4, colour="grey40") +
  geom_text(aes(label=sig, vjust=ifelse(delta>=0,-0.3,1.2)), size=2.8, fontface="bold") +
  scale_fill_manual(values=c(down=col_ctl, up=col_tbi), guide="none") +
  labs(title="D  Neuronal adhesion module, TBI vs control",
       subtitle="severity gradient; * = per-sample pseudobulk t-test; human shows strong collapse",
       x=NULL, y="Module score delta (TBI - control)") + pub

fig <- (pA | pB) / (pE | pC) / (pD | plot_spacer()) + plot_annotation(
  title="Computational support for the astrocyte-neuron adhesion axis (NRXN-led, cross-species)",
  theme=theme(plot.title=element_text(face="bold",size=11)))
ggsave(file.path(OUT,"FigS_computational_adhesion_axis.pdf"), fig, width=190, height=222, units="mm")
ggsave(file.path(OUT,"FigS_computational_adhesion_axis.png"), fig, width=190, height=222, units="mm", dpi=300)
cat("saved FigS_computational_adhesion_axis\n")
