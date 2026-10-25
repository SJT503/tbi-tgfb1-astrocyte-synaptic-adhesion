# JNI 手稿 — 发表级真实数据图统一重出（1F/6E + 补充图 S1/S3/S4/S5/S7）
# 统一发表主题、配色、字号；输出矢量 PDF + 300dpi PNG。数据=本项目真实分析输出。
suppressPackageStartupMessages({ library(ggplot2); library(grid); library(gridExtra) })

BASE <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project"
FH <- file.path(BASE,"wet_lab_data/results/figures/Fig_Human_TGFb")
FP <- file.path(BASE,"wet_lab_data/results/figures/Fig_PANoptosis")
FG <- file.path(BASE,"wet_lab_data/results/figures/Fig_Homeostatic_Deepdive")
FK <- file.path(BASE,"wet_lab_data/results/figures/Fig_Kir41_Gate")
FT <- file.path(BASE,"wet_lab_data/results/figures/Fig_TimeSeries")
FC <- file.path(BASE,"wet_lab_data/results/figures/Fig_CellChat")
PAN <- file.path(BASE,"JNI/figures/panels"); dir.create(PAN, showWarnings=FALSE, recursive=TRUE)
SUP <- file.path(BASE,"JNI/supplementary"); dir.create(SUP, showWarnings=FALSE, recursive=TRUE)

# ── 统一发表级主题 ──
pub <- theme_classic(base_size=9) +
  theme(text=element_text(family="sans", colour="black"),
        plot.title=element_text(face="bold", size=10, hjust=0, margin=margin(b=2)),
        plot.subtitle=element_text(size=7.4, colour="grey35", margin=margin(b=4)),
        axis.text=element_text(colour="black", size=8),
        axis.title=element_text(size=8.6, colour="black"),
        axis.line=element_line(colour="black", linewidth=0.35),
        axis.ticks=element_line(colour="black", linewidth=0.35),
        legend.position="top", legend.title=element_blank(),
        legend.text=element_text(size=7.6), legend.key.size=unit(3.2,"mm"),
        legend.margin=margin(b=-2),
        strip.background=element_rect(fill="#eef1f2", colour=NA),
        strip.text=element_text(face="bold", size=7.6, colour="#16324f"),
        plot.margin=margin(5,7,4,5))
col_ctl <- "#4DAECD"; col_tbi <- "#E07B54"
col_apop <- "#C0504D"; col_necr <- "#8B6BB1"
star <- function(p) ifelse(is.na(p),"",ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*",""))))
save_fig <- function(p, stem, w, h, dir=PAN){
  ggsave(file.path(dir,paste0(stem,".pdf")), p, width=w, height=h, units="mm")
  ggsave(file.path(dir,paste0(stem,".png")), p, width=w, height=h, units="mm", dpi=300)
}

# ════════ Fig 1F — 人源 TGFB1 各细胞类型 ════════
d <- read.csv(file.path(FH,"human_tgfb_byCelltype.csv"))
d1 <- d[d$gene=="TGFB1",]; d1$group <- factor(d1$group, levels=c("CTL","TBI"))
chg <- read.csv(file.path(FH,"human_TGFB1_change_byCelltype.csv")); chg$st <- star(chg$p)
ym <- tapply(d1$mean,d1$celltype,max)
sd1 <- data.frame(celltype=names(ym), y=as.numeric(ym)*1.13,
                  st=chg$st[match(names(ym),chg$celltype)])
p1F <- ggplot(d1, aes(celltype, mean, fill=group)) +
  geom_col(position=position_dodge(0.7), width=0.6, colour="black", linewidth=0.25) +
  geom_text(data=sd1, aes(celltype,y,label=st), inherit.aes=FALSE, size=3.6, fontface="bold") +
  scale_fill_manual(values=c(CTL=col_ctl,TBI=col_tbi)) +
  scale_y_continuous(expand=expansion(mult=c(0,.17))) +
  labs(title="Human snRNA-seq: TGFB1 by cell type",
       subtitle="Injured (TBI) vs control cortex, GSE209552", x=NULL, y="Mean expression") + pub
save_fig(p1F,"panel_1F_human_TGFB1", 84, 60)

# ════════ Fig 6E — 人源神经元死亡标志物 ════════
pa <- read.csv(file.path(FP,"panoptosis_external_neuron_panel.csv"))
h <- pa[pa$dataset=="GSE209552",]
sel <- c("Bax","Bcl2","Casp8","Ripk1","Mlkl")
arm <- c(Bax="Apoptosis",Bcl2="Apoptosis",Casp8="Apoptosis",Ripk1="Necroptosis",Mlkl="Necroptosis")
h <- h[h$gene %in% sel,]; h$arm <- arm[h$gene]; h$st <- star(h$p)
h$gene <- factor(h$gene, levels=sel)
h$lab <- ifelse(h$gene=="Bcl2","BCL2", toupper(as.character(h$gene)))
h$lab <- factor(h$lab, levels=h$lab[order(h$gene)])
p6E <- ggplot(h, aes(lab, pct, fill=arm)) +
  geom_col(width=0.6, colour="black", linewidth=0.25) +
  geom_hline(yintercept=0, linewidth=0.3) +
  geom_text(aes(label=st, vjust=ifelse(pct>=0,-0.25,1.15)), size=3.6, fontface="bold") +
  scale_fill_manual(values=c(Apoptosis=col_apop,Necroptosis=col_necr)) +
  scale_y_continuous(expand=expansion(mult=c(.1,.17))) +
  labs(title="Human snRNA-seq: neuronal death markers",
       subtitle="% change, injured (TBI) vs control cortex, GSE209552",
       x=NULL, y="% change (TBI vs control)") +
  pub + theme(axis.text.x=element_text(angle=25,hjust=1))
save_fig(p6E,"panel_6E_human_death", 88, 62)

# ════════ S1 — 单细胞 TGFβRI mRNA 跨细胞类型（5 数据集）════════
m <- read.csv(file.path(FH,"tgfbr1_byCelltype_external_mouse.csv"))
m <- data.frame(dataset=m$dataset, celltype=m$celltype, v=m$Tgfbr1_mean)
hu <- read.csv(file.path(FH,"human_TGFBR_change_byCelltype.csv"))
hu <- hu[hu$gene=="TGFBR1",]; hu <- data.frame(dataset="GSE209552 (human)", celltype=hu$celltype, v=hu$TBI)
sf <- read.csv(file.path(FC,"tgfbr_expression_by_celltype.csv")); names(sf)[1] <- "celltype"
nm <- c(Excitatory="Neuron",Inhibitory="Neuron",Astrocytes="Astro",Microglia="Micro",
        Oligodendrocytes="Oligo",OPCs="OPC",Pericyte="Endo")
sf2 <- data.frame(dataset="Self atlas",
        celltype=ifelse(sf$celltype %in% names(nm), nm[sf$celltype], sf$celltype), v=sf$Tgfbr1)
sf2 <- aggregate(v~dataset+celltype, sf2, max)
s1 <- rbind(m,hu,sf2)
s1$celltype <- factor(s1$celltype, levels=c("Micro","Endo","OPC","Oligo","Astro","Neuron"))
s1 <- s1[!is.na(s1$celltype),]
ctcol <- c(Micro="#C0504D",Endo="#9BBB59",OPC="#8064A2",Oligo="#4BACC6",Astro="#4DAECD",Neuron="#9aa6ab")
pS1 <- ggplot(s1, aes(celltype, v, fill=celltype)) +
  geom_col(width=0.74, colour="black", linewidth=0.22) +
  facet_wrap(~dataset, scales="free_y", nrow=1) +
  scale_fill_manual(values=ctcol) +
  scale_y_continuous(expand=expansion(mult=c(0,.08))) +
  labs(title="Supplementary Figure S1  |  TGFβRI mRNA across cell types in five single-cell datasets",
       subtitle="TGFβRI mRNA is microglia-high, whereas receptor protein localises to astrocytes (Fig. 2); mRNA does not equal protein localisation (see Discussion)",
       x=NULL, y="Mean expression") +
  pub + theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,size=7))
save_fig(pS1,"FigS1_TGFbRI_singlecell_mRNA", 250, 78, dir=SUP)

# ════════ S3 — PANoptosis 门控 ════════
agg <- aggregate(pct_det_TBI~gene+category, pa, function(v) mean(v,na.rm=TRUE))
agg$category <- factor(agg$category, levels=c("凋亡","坏死性凋亡","焦亡"),
                       labels=c("Apoptosis","Necroptosis","Pyroptosis"))
agg <- agg[order(agg$category,-agg$pct_det_TBI),]; agg$gene <- factor(agg$gene, levels=agg$gene)
pS3 <- ggplot(agg, aes(gene, pct_det_TBI, fill=category)) +
  geom_col(width=0.74, colour="black", linewidth=0.22) +
  geom_hline(yintercept=1, linetype="dashed", colour="#c0392b", linewidth=0.4) +
  facet_grid(~category, scales="free_x", space="free_x") +
  scale_fill_manual(values=c(Apoptosis=col_apop,Necroptosis=col_necr,Pyroptosis="#E0A030")) +
  scale_y_continuous(expand=expansion(mult=c(0,.08))) +
  labs(title="Supplementary Figure S3  |  Cell-death machinery detection in TBI neurons (single-cell)",
       subtitle="Mean of four datasets; necroptosis (Ripk3/Mlkl) and pyroptosis (Casp1/Gsdmd) executioner transcripts are near-undetectable (dashed line = 1%)",
       x=NULL, y="Detection rate in TBI neurons (%)") +
  pub + theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,size=7))
save_fig(pS3,"FigS3_PANoptosis_gating", 210, 80, dir=SUP)

# ════════ S4 — 谷氨酸/GS 面板（4 数据集）════════
gl <- read.csv(file.path(FG,"glutamate_panel_ALL_datasets.csv"))
gl <- gl[!gl$dataset %in% c("GSE292715","GSE230253","GSE226072"),]
gld <- rbind(
  data.frame(dataset=gl$dataset, gene="GLUL (GS)",      pct=gl$GLUL_pct,  p=gl$GLUL_p),
  data.frame(dataset=gl$dataset, gene="SLC1A2 (GLT-1)", pct=gl$SLC1A2_pct, p=gl$SLC1A2_p),
  data.frame(dataset=gl$dataset, gene="SLC1A3 (GLAST)", pct=gl$SLC1A3_pct, p=gl$SLC1A3_p))
gld$sig <- ifelse(gld$p<0.05,"P < 0.05","ns")
pS4 <- ggplot(gld, aes(dataset, pct, fill=sig)) +
  geom_col(width=0.66, colour="black", linewidth=0.22) +
  geom_hline(yintercept=0, linewidth=0.3) +
  facet_wrap(~gene, nrow=1) +
  scale_fill_manual(values=c("P < 0.05"=col_ctl,"ns"="grey78")) +
  labs(title="Supplementary Figure S4  |  Astrocyte glutamate-handling genes after TBI",
       subtitle="% change (TBI vs control), four datasets; GLUL/GS loss is robust, GLT-1/GLAST largely unchanged in mouse",
       x=NULL, y="% change (TBI vs control)") +
  pub + theme(axis.text.x=element_text(angle=45,hjust=1,size=7))
save_fig(pS4,"FigS4_glutamate_GS_panel", 220, 80, dir=SUP)

# ════════ S5 — 时间序列神经元 DEG ════════
FILES <- list(c("DEG_GSE271769_Neuron_D1.csv","GSE271769","D1"),
  c("DEG_GSE271769_Neuron_D7.csv","GSE271769","D7"),c("DEG_GSE271769_Neuron_D28.csv","GSE271769","D28"),
  c("DEG_GSE277487_Neuron_t6h.csv","GSE277487","6h"),c("DEG_GSE277487_Neuron_t2d.csv","GSE277487","2d"),
  c("DEG_GSE277487_Neuron_t4d.csv","GSE277487","4d"),
  c("DEG_GSE180862Hip_Neuron_24hrs.csv","GSE180862","24h"),
  c("DEG_GSE180862Hip_Neuron_7days.csv","GSE180862","7d"),
  c("DEG_GSE101901Hip_Neuron.csv","GSE101901","3d"))
deg <- do.call(rbind, lapply(FILES, function(x){
  fp <- file.path(FT,x[1]); if(!file.exists(fp)) return(NULL); d <- read.csv(fp)
  data.frame(dataset=x[2], tp=x[3],
    n_up=sum(d$padj<0.05 & d$log2FoldChange>1,na.rm=TRUE),
    n_down=sum(d$padj<0.05 & d$log2FoldChange< -1,na.rm=TRUE)) }))
degl <- rbind(data.frame(dataset=deg$dataset,tp=deg$tp,dir="Up-regulated",n=deg$n_up),
              data.frame(dataset=deg$dataset,tp=deg$tp,dir="Down-regulated",n=deg$n_down))
degl$tp <- factor(degl$tp, levels=c("6h","24h","2d","3d","D1","4d","7d","D7","D28"))
pS5 <- ggplot(degl, aes(tp,n,fill=dir)) +
  geom_col(width=0.66, colour="black", linewidth=0.22) +
  facet_wrap(~dataset, scales="free_x", nrow=1) +
  scale_fill_manual(values=c("Up-regulated"=col_tbi,"Down-regulated"=col_ctl)) +
  scale_y_continuous(expand=expansion(mult=c(0,.08))) +
  labs(title="Supplementary Figure S5  |  Neuronal transcriptional dysregulation over time after TBI",
       subtitle="Differentially expressed genes (pseudobulk DESeq2; |log2FC| > 1, adjusted P < 0.05) persist across post-injury time points",
       x=NULL, y="No. of DEGs") +
  pub + theme(strip.text=element_text(size=7.6))
save_fig(pS5,"FigS5_timeseries_neuron_DEG", 220, 76, dir=SUP)

# ════════ S7 — Kir4.1 阴性门禁 ════════
k <- read.csv(file.path(FK,"kir41_crossdataset_gate.csv"))
pa7 <- ggplot(k, aes(Aqp4_pct, Kir41_pct)) +
  geom_abline(slope=1,intercept=0,linetype="dashed",colour="grey60",linewidth=0.4) +
  geom_point(size=2.8, colour=col_necr) +
  geom_text(aes(label=dataset), size=2.3, vjust=-1.1, colour="grey35") +
  labs(title="A", subtitle="Kir4.1 tracks Aqp4: no Kir4.1-specific regulation",
       x="Aqp4 % change (TBI vs control)", y="Kir4.1 % change") + pub
pb7 <- ggplot(k, aes(reorder(dataset,Smad_Kir_r), Smad_Kir_r)) +
  geom_col(width=0.62, fill=col_necr, colour="black", linewidth=0.22) +
  geom_hline(yintercept=0,linewidth=0.3) + coord_flip() + ylim(-0.16,0.16) +
  labs(title="B", subtitle="Smad-Kir4.1 correlation near zero",
       x=NULL, y="Pearson r") + pub
S7 <- arrangeGrob(pa7,pb7,ncol=2,
  top=textGrob("Supplementary Figure S7  |  Kir4.1 fails cross-dataset gating (negative result)",
  x=0.02, hjust=0, gp=gpar(fontsize=10, fontface="bold", fontfamily="sans")))
ggsave(file.path(SUP,"FigS7_Kir41_negative.pdf"), S7, width=200, height=82, units="mm")
ggsave(file.path(SUP,"FigS7_Kir41_negative.png"), S7, width=200, height=82, units="mm", dpi=300)

cat("发表级真实数据图重出完成：1F/6E ->", PAN, " ; S1/S3/S4/S5/S7 ->", SUP, "\n")
