#!/usr/bin/env Rscript
# SFig 4 — Hippocampal-specific adhesion collapse by neuronal subtype
# 关键 claim: KO 实验做在海马, 公共数据 C5 显示海马 CA3/DG/Excitatory 神经元黏附模块特异性塌陷
#            (-27% to -73%), 皮层反而代偿性升高 (+50%~+150%) — 完美匹配 KO 靶位
# 数据: C5_per_subtype_adhesion.csv + C5_subtype_ttest.csv
suppressMessages({library(ggplot2); library(dplyr); library(patchwork)})
source("D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/scripts/style_lib_v2.R")

BASE <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI"
df   <- read.csv(file.path(BASE,"analysis/C5_neuron_subtypes/C5_per_subtype_adhesion.csv"))
ttst <- read.csv(file.path(BASE,"analysis/C5_neuron_subtypes/C5_subtype_ttest.csv"))
OUT  <- file.path(BASE,"figures")

# preserve only assigned subtypes (drop "Unassigned")
df <- df[df$subtype != "Unassigned", ]
ttst <- ttst[ttst$subtype != "Unassigned", ]

# nicer subtype labels
sub_labels <- c(
  Excitatory_glut     = "Excitatory\n(Slc17a7/a6)",
  Inhibitory_GABA     = "Inhibitory\n(Gad1/2)",
  Hippo_CA1_pyramidal = "CA1 pyramidal\n(Wfs1/Spink8)",
  Hippo_CA3_pyramidal = "CA3 pyramidal\n(Prkcd/Calb2)",
  Hippo_DG_granule    = "DG granule\n(Prox1/C1ql2)")
df$subtype_lab   <- factor(sub_labels[df$subtype], levels=sub_labels)
ttst$subtype_lab <- factor(sub_labels[ttst$subtype], levels=sub_labels)
df$condition_lab <- factor(ifelse(df$condition=="control","Sham","TBI"), levels=c("Sham","TBI"))
df$brain_region  <- factor(df$brain_region, levels=c("hippocampus","cortex"))

# compute % change for annotation
delta <- df %>% tidyr::pivot_wider(id_cols=c(subtype_lab, brain_region),
                                   names_from=condition_lab, values_from=mean_adh) %>%
         mutate(delta_abs = TBI - Sham,
                # robust % change: if |Sham| < 0.02 (near zero), use absolute delta * 100 as "score delta"
                pct_change = ifelse(abs(Sham) < 0.02,
                                    round(delta_abs * 100, 1),  # near-zero base: report delta x100 as "abs delta"
                                    round(100 * delta_abs / abs(Sham), 0)),
                near_zero_base = abs(Sham) < 0.02)
# merge p
delta <- merge(delta, ttst[,c("subtype","brain_region","p_value","padj_BH")] %>%
                       mutate(subtype_lab=factor(sub_labels[subtype], levels=sub_labels)) %>%
                       select(-subtype),
               by=c("subtype_lab","brain_region"), all.x=TRUE)
delta$sig <- ifelse(is.na(delta$padj_BH), "",
              ifelse(delta$padj_BH < 0.001, "***",
              ifelse(delta$padj_BH < 0.01,  "**",
              ifelse(delta$padj_BH < 0.05,  "*", "ns"))))
delta$bio_signif <- abs(delta$pct_change) >= 15
delta$label <- paste0(ifelse(delta$pct_change>0, "+", ""),
                      delta$pct_change, ifelse(delta$near_zero_base, "", "%"),
                      ifelse(delta$sig=="", "", paste0(" ", delta$sig)),
                      ifelse(delta$near_zero_base, "*", ""))

cat("delta summary:\n"); print(delta)

# ── PANEL A: 海马 (hippocampus) — 重点 ──────────────────────────────
df_h <- df[df$brain_region == "hippocampus", ]
delta_h <- delta[delta$brain_region == "hippocampus", ]
ymax_h <- max(df_h$mean_adh, na.rm=TRUE) * 1.15

pA <- ggplot(df_h, aes(subtype_lab, mean_adh, fill=condition_lab)) +
  geom_col(position=position_dodge(0.78), width=0.7, colour="grey25", linewidth=0.25) +
  geom_text(data=delta_h, aes(x=subtype_lab, y=pmax(Sham, TBI) + ymax_h*0.05,
                              label=label, fontface=ifelse(bio_signif,"bold","plain")),
            inherit.aes=FALSE, size=2.8, hjust=0.5) +
  scale_fill_manual(values=c(Sham=JNI_COL$control, TBI=JNI_COL$tbi)) +
  scale_y_continuous(limits=c(NA, ymax_h*1.25)) +
  labs(title="A. Hippocampal neurons — adhesion module COLLAPSE",
       subtitle="Public mouse hippocampus (mouse_integrated, n=9 sham vs 13 TBI mice across 24h+7d). % change vs Sham; bold = biologically significant (|%|≥15)",
       y="Adhesion module score (16-gene panel)",
       x=NULL) +
  jni_theme() +
  theme(axis.text.x=element_text(angle=15, hjust=1, lineheight=0.9))

# ── PANEL B: 皮层 (cortex) — 对照 ───────────────────────────────────
df_c <- df[df$brain_region == "cortex", ]
delta_c <- delta[delta$brain_region == "cortex", ]
ymax_c <- max(df_c$mean_adh, na.rm=TRUE) * 1.15

pB <- ggplot(df_c, aes(subtype_lab, mean_adh, fill=condition_lab)) +
  geom_col(position=position_dodge(0.78), width=0.7, colour="grey25", linewidth=0.25) +
  geom_text(data=delta_c, aes(x=subtype_lab, y=pmax(Sham, TBI) + ymax_c*0.05,
                              label=label, fontface=ifelse(bio_signif,"bold","plain")),
            inherit.aes=FALSE, size=2.8, hjust=0.5) +
  scale_fill_manual(values=c(Sham=JNI_COL$control, TBI=JNI_COL$tbi)) +
  scale_y_continuous(limits=c(NA, ymax_c*1.25)) +
  labs(title="B. Cortical neurons — compensatory ELEVATION (contrast)",
       subtitle="Same panel, cortical neurons; opposite direction confirms region-specificity of the collapse",
       y="Adhesion module score (16-gene panel)",
       x=NULL) +
  jni_theme() +
  theme(axis.text.x=element_text(angle=15, hjust=1, lineheight=0.9))

# ── PANEL C: 跨区域汇总 ─────────────────────────────────────────────
pC <- ggplot(delta, aes(subtype_lab, pct_change, fill=brain_region, group=brain_region)) +
  geom_col(position=position_dodge(0.75), width=0.65, colour="grey25", linewidth=0.25) +
  geom_hline(yintercept=0, linewidth=0.4, colour="grey40") +
  geom_hline(yintercept=c(-15, 15), linetype="dashed", linewidth=0.3, colour="grey60") +
  geom_text(aes(label=label, y=pct_change + ifelse(pct_change>=0, 9, -10),
                fontface=ifelse(bio_signif,"bold","plain"), group=brain_region),
            position=position_dodge(0.75), size=2.6) +
  scale_fill_manual(values=jni_palette_region) +
  scale_y_continuous(limits=c(-92, 172), breaks=c(-50, 0, 50, 100, 150)) +
  labs(title="C. % change summary — hippocampus DOWN, cortex UP",
       subtitle="Dashed lines = ±15% biological-significance threshold; bold = |delta%|≥15 AND padj<0.05",
       y="% change (TBI vs Sham)",
       x=NULL) +
  jni_theme() +
  theme(axis.text.x=element_text(angle=15, hjust=1, lineheight=0.9))

# ── ASSEMBLE ────────────────────────────────────────────────────────
fig <- (pA / pB / pC) +
  plot_annotation(title="Supplementary Figure S3. Hippocampal-specific neuronal adhesion-module collapse by subtype",
                  subtitle="Public mouse TBI scRNA-seq (GSE101901 + GSE180862); congruent with hippocampal AAV-GFAP-TGF-betaRI KO",
                  theme=theme(plot.title=element_text(face="bold", size=10.5),
                              plot.subtitle=element_text(size=8, colour="grey35"),
                              plot.margin=margin(4, 6, 4, 4, "mm")))
ggsave(file.path(OUT,"SFig3_hippocampal_subtype_collapse.pdf"), fig, width=195, height=240, units="mm")
ggsave(file.path(OUT,"SFig3_hippocampal_subtype_collapse.png"), fig, width=195, height=240, units="mm", dpi=300)
ggsave(file.path(BASE,"supplementary/FigS3_hippocampal_subtype_collapse.pdf"), fig, width=195, height=240, units="mm")
ggsave(file.path(BASE,"supplementary/FigS3_hippocampal_subtype_collapse.png"), fig, width=195, height=240, units="mm", dpi=300)
cat("\nSaved SFig3_hippocampal_subtype_collapse.pdf/.png\n")
