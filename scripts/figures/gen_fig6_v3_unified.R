#!/usr/bin/env Rscript
# Fig 6 (v3, R2 升级) — 用 style_lib_v2.R 统一风格 + 加 effect size 注释
# 仍是 4-panel 黏附塌陷主图: A 小鼠 CellChat / B 人源 CellChat / C 自测三组 / D 神经元黏附跨严重度
# v3 改进:
#   - 调色板从 hard-coded 切换到 jni_palette_3 / jni_palette_2 / JNI_COL
#   - 主题用 jni_theme() 统一字体/字号/grid
#   - panel 标题用 capital letter + bold
#   - effect size 注释（% change）直接打在柱顶
suppressMessages({library(ggplot2); library(patchwork); library(dplyr)})
source("D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/scripts/style_lib_v2.R")

BASE <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project"
S5 <- file.path(BASE,"JNI/analysis/step5_out")
S6 <- file.path(BASE,"JNI/analysis/step6_out")
CO <- file.path(BASE,"JNI/analysis/cand_out")
OUT<- file.path(BASE,"JNI/figures")

# ── PANEL A: 小鼠公共 CellChat astrocyte→neuron 黏附通路 ────────────
b <- read.csv(file.path(S5,"5b_astro_to_neuron_pathways.csv"))
adh <- c("NRXN","NCAM","NEGR","NGL","CNTN","CADM")
b <- b[b$pathway_name %in% adh,]
b$pathway_name <- factor(b$pathway_name, levels=adh)
b$dir <- ifelse(b$log2FC<0,"down","up")
b$pct <- round(100 * (2^b$log2FC - 1), 0)  # log2FC → % change
b$label <- ifelse(b$pct == 0, "0%",
            paste0(ifelse(b$pct>0,"+",""), b$pct, "%"))
pA <- ggplot(b, aes(pathway_name, log2FC, fill=dir)) +
  geom_col(width=0.62, colour="grey25", linewidth=0.25) +
  geom_hline(yintercept=0, linewidth=0.4, colour="grey40") +
  geom_text(aes(label=label, y=log2FC + ifelse(log2FC>=0, 0.055, -0.065)),
            size=2.6, fontface="plain") +
  scale_fill_manual(values=c(down=JNI_COL$down, up=JNI_COL$up), guide="none") +
  scale_y_continuous(limits=c(-1.08, 0.18), breaks=seq(-1.0, 0.0, by=0.25)) +
  labs(title="A  Mouse: astrocyte-to-neuron adhesion (CellChat)",
       subtitle="Public mouse TBI cortex; NRXN/NEGR/CNTN/NGL collapse after injury",
       x=NULL, y="log2(TBI / control)") +
  jni_theme(base=9)

# ── PANEL B: 人源 CellChat (cand2b strict) ──────────────────────────
e <- read.csv(file.path(CO,"cand2b_human_astro2neuron_pathways.csv"))
e <- e[e$adhesion==TRUE & !is.na(e$pct_change) & is.finite(e$pct_change),]
e$pathway <- factor(e$pathway, levels=e$pathway[order(e$pct_change)])
e$dir <- ifelse(e$pct_change<0,"down","up")
e$label <- paste0(ifelse(e$pct_change>0,"+",""), round(e$pct_change, 0), "%")
pB <- ggplot(e, aes(pathway, pct_change, fill=dir)) +
  geom_col(width=0.62, colour="grey25", linewidth=0.25) +
  geom_hline(yintercept=0, linewidth=0.4, colour="grey40") +
  geom_text(aes(label=label, y=pct_change + ifelse(pct_change>=0, 6, -6.5)),
            size=2.6, fontface="plain") +
  scale_fill_manual(values=c(down=JNI_COL$down, up=JNI_COL$up), guide="none") +
  scale_y_continuous(limits=c(-114, 48), breaks=c(-100, -50, 0, 40)) +
  labs(title="B  Human: astrocyte-to-neuron adhesion (CellChat)",
       subtitle="GSE209552 TBI cortex; NRXN-led collapse, NCAM discordant",
       x=NULL, y="% change (TBI vs control)") +
  jni_theme(base=9)

# ── PANEL C: 自测三组 CellChat ─────────────────────────────────────
cc <- read.csv(file.path(BASE,"wet_lab_data/results/cellchat/Ast_to_DG_pathway_strength_3groups.csv"))
cc <- cc[,c("Pathway","Sham","CCI_Scr","KO")]
cl <- reshape(cc, varying=c("Sham","CCI_Scr","KO"), v.names="strength",
              timevar="group", times=c("Sham","CCI_Scr","KO"), direction="long")
cl$group <- factor(cl$group, levels=c("Sham","CCI_Scr","KO"))
cl$Pathway <- factor(cl$Pathway, levels=cc$Pathway[order(-cc$Sham)])
pC <- ggplot(cl, aes(Pathway, strength, fill=group)) +
  geom_col(position=position_dodge(0.76), width=0.7,
           colour="grey25", linewidth=0.25) +
  scale_fill_manual(values=c(Sham=JNI_COL$control,
                              CCI_Scr=JNI_COL$tbi, KO=JNI_COL$ko)) +
  labs(title="C  Self snRNA three groups (descriptive, n=1/group)",
       subtitle="CCI lowers, astrocytic TGFBRI-KO restores; descriptive only",
       x=NULL, y="Communication strength") +
  jni_theme(base=9, legend="top")

# ── PANEL D: 神经元黏附模块跨严重度 ─────────────────────────────────
d <- read.csv(file.path(S6,"step6b_severity_gradient.csv"))
d <- d[d$module=="Adhesion",]
lab <- c(GSE290150="Moderate",GSE271769="Severe-1",
         GSE277487="Severe-2",GSE209552="Human")
d$dataset <- factor(lab[d$dataset], levels=lab)
d$dir <- ifelse(d$delta<0,"down","up")
d$sig <- ifelse(d$pseudobulk_p<0.001,"***",
          ifelse(d$pseudobulk_p<0.01,"**",
          ifelse(d$pseudobulk_p<0.05,"*","ns")))
d$pct <- round(100 * d$delta / abs(pmax(0.01, abs(d$delta) + 0.01)), 0)  # approximate % change
# better: use relative module change; here delta is already in module-score units
d$label <- paste0(ifelse(d$delta>=0,"+",""), round(d$delta, 3), " ", d$sig)
pD <- ggplot(d, aes(dataset, delta, fill=dir)) +
  geom_col(width=0.6, colour="grey25", linewidth=0.25) +
  geom_hline(yintercept=0, linewidth=0.4, colour="grey40") +
  geom_text(aes(label=label,
                y=delta + ifelse(delta>=0, 0.028, -0.032)),
            size=2.6, fontface="bold") +
  scale_fill_manual(values=c(down=JNI_COL$down, up=JNI_COL$up), guide="none") +
  scale_y_continuous(limits=c(-0.49, 0.155), breaks=c(-0.4, -0.3, -0.2, -0.1, 0.0, 0.1)) +
  labs(title="D  Neuron adhesion module across severity (cortex)",
       subtitle="Cortical neurons across severity; human shows strongest collapse",
       x=NULL, y="Module score delta (TBI - control)") +
  jni_theme(base=9)

# ── ASSEMBLE ────────────────────────────────────────────────────────
fig <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title="Figure 6. Astrocytic TGF-betaRI reprograms the astrocyte-neuron synaptic-adhesion programme",
    subtitle="NRXN-led adhesion collapse, conserved across mouse and human TBI cortex; restored by astrocyte-specific TGF-betaRI knockout",
    theme=theme(plot.title=element_text(face="bold", size=10.5),
                plot.subtitle=element_text(size=8.5, colour="grey35"),
                plot.margin=margin(4, 6, 4, 4, "mm")))

ggsave(file.path(OUT,"Fig6_v3_unified.pdf"), fig, width=195, height=152, units="mm")
ggsave(file.path(OUT,"Fig6_v3_unified.png"), fig, width=195, height=152, units="mm", dpi=300)
ggsave(file.path(BASE,"JNI/课题最新框架和待补实验/Fig6.pdf"), fig, width=195, height=152, units="mm")
cat("saved Fig6_v3_unified (style_lib_v2.R)\n")
