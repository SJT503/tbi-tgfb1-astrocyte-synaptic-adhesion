# JNI 手稿 — 补充图 S1/S3/S4/S5/S7 生成（真实数据，来自已有分析脚本 CSV 输出）
suppressPackageStartupMessages({ library(ggplot2) })
BASE <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project"
FH   <- file.path(BASE,"wet_lab_data/results/figures/Fig_Human_TGFb")
FP   <- file.path(BASE,"wet_lab_data/results/figures/Fig_PANoptosis")
FG   <- file.path(BASE,"wet_lab_data/results/figures/Fig_Homeostatic_Deepdive")
FK   <- file.path(BASE,"wet_lab_data/results/figures/Fig_Kir41_Gate")
FT   <- file.path(BASE,"wet_lab_data/results/figures/Fig_TimeSeries")
FC   <- file.path(BASE,"wet_lab_data/results/figures/Fig_CellChat")
OUT  <- file.path(BASE,"JNI/supplementary")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
thm <- theme_classic(base_size=10) +
  theme(text=element_text(family="sans"), plot.title=element_text(face="bold",size=10.5),
        plot.subtitle=element_text(size=8,color="#666"),
        axis.text=element_text(color="black"), legend.position="top",
        legend.title=element_blank(), legend.key.size=unit(3.5,"mm"))

# ───────── S1 — 单细胞 TGFβRI mRNA 跨细胞类型分布（5 数据集）─────────
m  <- read.csv(file.path(FH,"tgfbr1_byCelltype_external_mouse.csv"))
m  <- data.frame(dataset=m$dataset, celltype=m$celltype, tgfbr1=m$Tgfbr1_mean)
hu <- read.csv(file.path(FH,"human_TGFBR_change_byCelltype.csv"))
hu <- hu[hu$gene=="TGFBR1",]; hu <- data.frame(dataset="GSE209552(human)",
        celltype=hu$celltype, tgfbr1=hu$TBI)
sf <- read.csv(file.path(FC,"tgfbr_expression_by_celltype.csv"))
names(sf)[1] <- "celltype"
sf2 <- data.frame(dataset="Self atlas", celltype=sf$celltype, tgfbr1=sf$Tgfbr1)
norm <- c(Excitatory="Neuron",Inhibitory="Neuron",Astrocytes="Astro",
          Microglia="Micro",Oligodendrocytes="Oligo",OPCs="OPC",Pericyte="Endo")
sf2$celltype <- ifelse(sf2$celltype %in% names(norm), norm[sf2$celltype], sf2$celltype)
sf2 <- aggregate(tgfbr1~dataset+celltype, sf2, max)
s1 <- rbind(m, hu, sf2)
s1$celltype <- factor(s1$celltype, levels=c("Micro","Endo","OPC","Oligo","Astro","Neuron"))
s1 <- s1[!is.na(s1$celltype),]
pS1 <- ggplot(s1, aes(celltype, tgfbr1, fill=celltype)) +
  geom_col(width=0.7) + facet_wrap(~dataset, scales="free_y", nrow=1) +
  scale_fill_manual(values=c(Micro="#C0504D",Endo="#9BBB59",OPC="#8064A2",
    Oligo="#4BACC6",Astro="#4DAECD",Neuron="#808080")) +
  labs(title="Supplementary Figure S1   TGFβRI (Tgfbr1/TGFBR1) mRNA across cell types — 5 single-cell datasets",
       subtitle="TGFbRI mRNA is microglia-high in single-cell data, whereas protein IF localizes to astrocytes (main Fig 2): mRNA does not equal protein localization (see Discussion)",
       x=NULL, y="Mean expression") +
  thm + theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,size=7.5),
        strip.text=element_text(size=7.6,face="bold"))
ggsave(file.path(OUT,"FigS1_TGFbRI_singlecell_mRNA.pdf"), pS1, width=11, height=3.4)

# ───────── S3 — PANoptosis 单细胞门控（神经元死亡机器检出率）─────────
pa <- read.csv(file.path(FP,"panoptosis_external_neuron_panel.csv"))
agg <- aggregate(pct_det_TBI~gene+category, pa, function(v) mean(v,na.rm=TRUE))
catorder <- c("凋亡","坏死性凋亡","焦亡")
agg$category <- factor(agg$category, levels=catorder,
  labels=c("Apoptosis","Necroptosis","Pyroptosis"))
agg <- agg[order(agg$category, -agg$pct_det_TBI),]
agg$gene <- factor(agg$gene, levels=agg$gene)
pS3 <- ggplot(agg, aes(gene, pct_det_TBI, fill=category)) +
  geom_col(width=0.72) +
  geom_hline(yintercept=1, linetype="dashed", color="#c0392b", linewidth=0.4) +
  facet_grid(~category, scales="free_x", space="free_x") +
  scale_fill_manual(values=c(Apoptosis="#C0504D",Necroptosis="#8B6BB1",Pyroptosis="#E0A030")) +
  labs(title="Supplementary Figure S3   PANoptosis machinery detection in TBI neurons (single-cell, mean of 4 datasets)",
       subtitle="Necroptosis (Ripk3/Mlkl) and pyroptosis (Casp1/Gsdmd) executioner genes are near-undetectable in neuron transcriptomes (dashed line = 1%): a transcriptomic limitation, not a negative result",
       x=NULL, y="Detection rate in TBI neurons (%)") +
  thm + theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,size=7.5),
        strip.text=element_text(size=8,face="bold"))
ggsave(file.path(OUT,"FigS3_PANoptosis_gating.pdf"), pS3, width=9, height=3.6)

# ───────── S4 — 谷氨酸/GS 面板（删除 3 个已舍弃数据集）─────────
gl <- read.csv(file.path(FG,"glutamate_panel_ALL_datasets.csv"))
drop <- c("GSE292715","GSE230253","GSE226072")
gl <- gl[!gl$dataset %in% drop,]
gld <- rbind(
  data.frame(dataset=gl$dataset, gene="GLUL (GS)",  pct=gl$GLUL_pct,  p=gl$GLUL_p),
  data.frame(dataset=gl$dataset, gene="SLC1A2 (GLT-1)", pct=gl$SLC1A2_pct, p=gl$SLC1A2_p),
  data.frame(dataset=gl$dataset, gene="SLC1A3 (GLAST)", pct=gl$SLC1A3_pct, p=gl$SLC1A3_p))
gld$sig <- ifelse(gld$p<0.05,"P<0.05","ns")
pS4 <- ggplot(gld, aes(dataset, pct, fill=sig)) +
  geom_col(width=0.66) + geom_hline(yintercept=0, linewidth=0.3) +
  facet_wrap(~gene, nrow=1) +
  scale_fill_manual(values=c("P<0.05"="#4DAECD","ns"="grey75")) +
  labs(title="Supplementary Figure S4   Astrocyte glutamate-handling genes — TBI vs control (% change), 4 datasets",
       subtitle="GLUL/GS loss is robust; GLT-1/GLAST largely unchanged in mouse: astrocyte reactive-state context only (3 discarded datasets removed)",
       x=NULL, y="% change (TBI vs CTL)") +
  thm + theme(axis.text.x=element_text(angle=45,hjust=1,size=7.5),
        strip.text=element_text(size=8,face="bold"))
ggsave(file.path(OUT,"FigS4_glutamate_GS_panel.pdf"), pS4, width=10, height=3.6)

# ───────── S5 — 时间序列神经元 DEG 累积 ─────────
FILES <- list(
  c("DEG_GSE271769_Neuron_D1.csv","GSE271769","D1"),
  c("DEG_GSE271769_Neuron_D7.csv","GSE271769","D7"),
  c("DEG_GSE271769_Neuron_D28.csv","GSE271769","D28"),
  c("DEG_GSE277487_Neuron_t6h.csv","GSE277487","6h"),
  c("DEG_GSE277487_Neuron_t2d.csv","GSE277487","2d"),
  c("DEG_GSE277487_Neuron_t4d.csv","GSE277487","4d"),
  c("DEG_GSE180862Hip_Neuron_24hrs.csv","GSE180862","24h"),
  c("DEG_GSE180862Hip_Neuron_7days.csv","GSE180862","7d"),
  c("DEG_GSE101901Hip_Neuron.csv","GSE101901","3d"))
deg <- do.call(rbind, lapply(FILES, function(x){
  fp <- file.path(FT,x[1]); if(!file.exists(fp)) return(NULL)
  d <- read.csv(fp)
  data.frame(dataset=x[2], tp=x[3],
    n_up  =sum(d$padj<0.05 & d$log2FoldChange> 1, na.rm=TRUE),
    n_down=sum(d$padj<0.05 & d$log2FoldChange< -1, na.rm=TRUE)) }))
degl <- rbind(data.frame(dataset=deg$dataset, tp=deg$tp, dir="Up",   n=deg$n_up),
              data.frame(dataset=deg$dataset, tp=deg$tp, dir="Down", n=deg$n_down))
degl$tp <- factor(degl$tp, levels=c("6h","24h","2d","3d","D1","4d","7d","D7","D28"))
pS5 <- ggplot(degl, aes(tp, n, fill=dir)) +
  geom_col(width=0.66) + facet_wrap(~dataset, scales="free_x", nrow=1) +
  scale_fill_manual(values=c(Up="#E07B54",Down="#4DAECD")) +
  labs(title="Supplementary Figure S5   Neuronal transcriptional dysregulation over time (pseudobulk DESeq2 DEGs, |log2FC|>1, padj<0.05)",
       subtitle="Neuronal DEG counts persist across post-TBI time points: chronic transcriptional dysregulation",
       x=NULL, y="No. of DEGs") +
  thm + theme(strip.text=element_text(size=7.6,face="bold"))
ggsave(file.path(OUT,"FigS5_timeseries_neuron_DEG.pdf"), pS5, width=10, height=3.4)

# ───────── S7 — Kir4.1 阴性门禁 ─────────
k <- read.csv(file.path(FK,"kir41_crossdataset_gate.csv"))
pS7a <- ggplot(k, aes(Aqp4_pct, Kir41_pct)) +
  geom_abline(slope=1,intercept=0,linetype="dashed",color="grey60") +
  geom_point(size=2.6, color="#8B6BB1") +
  geom_text(aes(label=dataset), size=2.3, vjust=-1, color="#555") +
  labs(title="A  Kir4.1 vs Aqp4 % change", subtitle="Kir4.1 tracks Aqp4: no Kir4.1-specific regulation",
       x="Aqp4 % change", y="Kir4.1 % change") + thm
pS7b <- ggplot(k, aes(reorder(dataset,Smad_Kir_r), Smad_Kir_r)) +
  geom_col(width=0.6, fill="#8B6BB1") + geom_hline(yintercept=0,linewidth=0.3) +
  coord_flip() + ylim(-0.15,0.15) +
  labs(title="B  Smad-Kir4.1 correlation", subtitle="Correlation near zero: Kir4.1 not downstream of TGF-b/Smad",
       x=NULL, y="Pearson r") + thm
pdf(file.path(OUT,"FigS7_Kir41_negative.pdf"), width=9, height=3.6)
gridExtra::grid.arrange(pS7a, pS7b, ncol=2,
  top=grid::textGrob("Supplementary Figure S7   Kir4.1 fails cross-dataset gating (negative result)",
  gp=grid::gpar(fontsize=11, fontface="bold")))
dev.off()

cat("补充图 S1/S3/S4/S5/S7 已生成 ->", OUT, "\n")
