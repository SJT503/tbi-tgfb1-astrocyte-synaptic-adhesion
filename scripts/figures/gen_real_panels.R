# JNI 手稿 — 真实数据 panel 生成（Fig 1F 人源 TGFB1；Fig 6E 人源双死亡标志物）
# 数据来源：本项目已有分析脚本输出 CSV（真实数据，非模拟）
suppressPackageStartupMessages({ library(ggplot2) })
BASE <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project"
FIGD <- file.path(BASE,"wet_lab_data/results/figures/Fig_Human_TGFb")
PANO <- file.path(BASE,"wet_lab_data/results/figures/Fig_PANoptosis")
OUT  <- file.path(BASE,"JNI/figures/panels")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

# JNI 配色
col_ctl <- "#4DAECD"; col_tbi <- "#E07B54"
thm <- theme_classic(base_size=10) +
  theme(text=element_text(family="sans"),
        plot.title=element_text(face="bold", size=10.5),
        axis.text=element_text(color="black"),
        legend.position="top", legend.title=element_blank(),
        legend.key.size=unit(3.5,"mm"))

# ───────── Fig 1F — 人源 snRNA TGFB1 各细胞类型 TBI vs CTL ─────────
d <- read.csv(file.path(FIGD,"human_tgfb_byCelltype.csv"))
d1 <- d[d$gene=="TGFB1",]
d1$group <- factor(d1$group, levels=c("CTL","TBI"))
chg <- read.csv(file.path(FIGD,"human_TGFB1_change_byCelltype.csv"))
chg$star <- ifelse(chg$p<0.001,"***",ifelse(chg$p<0.01,"**",ifelse(chg$p<0.05,"*","")))
ymax <- tapply(d1$mean, d1$celltype, max)
star_df <- data.frame(celltype=names(ymax), y=as.numeric(ymax)*1.12,
                      star=chg$star[match(names(ymax), chg$celltype)])

p1F <- ggplot(d1, aes(celltype, mean, fill=group)) +
  geom_col(position=position_dodge(0.72), width=0.64) +
  geom_text(data=star_df, aes(celltype, y, label=star), inherit.aes=FALSE,
            size=4.2, fontface="bold") +
  scale_fill_manual(values=c(CTL=col_ctl, TBI=col_tbi)) +
  labs(title="F   Human snRNA: TGFB1 by cell type (TBI vs control)",
       x=NULL, y="Mean expression") +
  scale_y_continuous(expand=expansion(mult=c(0,.16))) + thm
ggsave(file.path(OUT,"panel_1F_human_TGFB1.png"), p1F, width=82, height=62,
       units="mm", dpi=300)
cat("Fig 1F 已出\n")

# ───────── Fig 6E — 人源 snRNA 神经元双死亡标志物（TBI vs CTL %变化）─────────
pa <- read.csv(file.path(PANO,"panoptosis_external_neuron_panel.csv"))
h <- pa[pa$dataset=="GSE209552",]
sel <- c("Bax","Bcl2","Casp8","Ripk1","Mlkl")          # 凋亡 + 坏死性凋亡 curated
arm <- c(Bax="Apoptosis",Bcl2="Apoptosis",Casp8="Apoptosis",
         Ripk1="Necroptosis",Mlkl="Necroptosis")
h <- h[h$gene %in% sel,]
h$arm <- arm[h$gene]
h$star <- ifelse(h$p<0.001,"***",ifelse(h$p<0.01,"**",ifelse(h$p<0.05,"*","")))
h$gene <- factor(h$gene, levels=sel)
h$lab <- ifelse(h$gene=="Bcl2","BCL2 (anti-apop)", toupper(as.character(h$gene)))
h$lab <- factor(h$lab, levels=h$lab[order(h$gene)])

p6E <- ggplot(h, aes(lab, pct, fill=arm)) +
  geom_col(width=0.62) +
  geom_hline(yintercept=0, linewidth=0.3) +
  geom_text(aes(label=star, vjust=ifelse(pct>=0,-0.2,1.1)), size=4, fontface="bold") +
  scale_fill_manual(values=c(Apoptosis="#C0504D", Necroptosis="#8B6BB1")) +
  labs(title="E   Human snRNA: neuronal death markers (TBI vs control)",
       x=NULL, y="% change (TBI vs CTL)") +
  scale_y_continuous(expand=expansion(mult=c(.08,.16))) +
  thm + theme(axis.text.x=element_text(angle=20, hjust=1, size=8.5))
ggsave(file.path(OUT,"panel_6E_human_death.png"), p6E, width=86, height=64,
       units="mm", dpi=300)
cat("Fig 6E 已出\n")
cat("\n真实数据 panel 完成 ->", OUT, "\n")
