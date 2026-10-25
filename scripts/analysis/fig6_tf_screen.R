#!/usr/bin/env Rscript
# 研究基础图6 —— 无偏 TF 筛选锁定 FOXG1 (JNI 风格, 真实筛选数据)
suppressMessages({library(ggplot2); library(cowplot)})
P <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/tf_screen_out"
df <- read.csv(file.path(P,"CROSS_3datasets_TFscreen.csv"), stringsAsFactors=FALSE)
cat("TF 总数:", nrow(df), " 候选②(n_lostact==3):", sum(df$n_lostact==3), "\n")
fx <- df[df$TF=="FOXG1",]
stopifnot(nrow(fx)==1)
cat("FOXG1: mean_lfc=",fx$mean_lfc," mean_rho=",fx$mean_rho,"\n")

navy<-"#1c3553"; teal<-"#2f8f88"; gold<-"#b8932f"
df$grp <- ifelse(df$TF=="FOXG1","FOXG1", ifelse(df$n_lostact==3,"cand","other"))

# ---- Panel A: 散点 ----
pA <- ggplot(df, aes(mean_lfc, mean_rho)) +
  geom_hline(yintercept=0, linetype=2, color="grey75") +
  geom_vline(xintercept=0, linetype=2, color="grey75") +
  geom_point(data=df[df$grp=="other",], color="grey82", size=1.5) +
  geom_point(data=df[df$grp=="cand",], color=teal, size=2.3, alpha=.9) +
  geom_point(data=df[df$grp=="FOXG1",], color=gold, size=4.6) +
  geom_text(data=df[df$grp=="FOXG1",], aes(label=TF), vjust=-1.25,
            color=navy, fontface="bold", size=4.4) +
  annotate("text", x=min(df$mean_lfc)*0.96, y=max(df$mean_rho)*0.97,
           label="injury-lost, adhesion-supporting TFs", hjust=0,
           color=teal, fontface="italic", size=3.5) +
  labs(x="Injury effect  (mean log2FC, TBI vs control)",
       y="Correlation with synaptic-adhesion program (mean Spearman ρ)",
       title="A   Unbiased TF screen across 3 independent datasets") +
  theme_cowplot(11) + theme(plot.title=element_text(size=11.5,color=navy))

# ---- Panel B/C: FOXG1 跨数据集 ----
ds <- c("GSE226207\n(mouse)","GSE101901/180862\n(mouse)","GSE209552\n(human)")
b <- data.frame(ds=factor(ds,levels=ds),
                lfc=c(fx$lfc_226207,fx$lfc_integ,fx$lfc_209552),
                rho=c(fx$rho_226207,fx$rho_integ,fx$rho_209552))
pB <- ggplot(b, aes(ds,lfc)) +
  geom_col(fill=navy, width=.62) +
  geom_hline(yintercept=0,color="grey40") +
  geom_text(aes(label=sprintf("%.2f",lfc)), vjust=1.35, color=navy, size=3.5, fontface="bold") +
  scale_y_continuous(expand=expansion(mult=c(.20,.06))) +
  labs(x=NULL, y="FOXG1 log2FC",
       title="B   FOXG1 consistently down-regulated after injury") +
  theme_cowplot(11) + theme(plot.title=element_text(size=11,color=navy),
       axis.text.x=element_text(size=8.3),
       plot.margin=margin(6,8,4,8))
pC <- ggplot(b, aes(ds,rho)) +
  geom_col(fill=teal, width=.62) +
  geom_text(aes(label=sprintf("%.2f",rho)), vjust=-0.5, color=navy, size=3.5, fontface="bold") +
  scale_y_continuous(expand=expansion(mult=c(.04,.20))) +
  labs(x=NULL, y="FOXG1–adhesion ρ",
       title="C   FOXG1 positively tracks the adhesion program") +
  theme_cowplot(11) + theme(plot.title=element_text(size=11,color=navy),
       axis.text.x=element_text(size=8.3),
       plot.margin=margin(6,8,4,8))

right <- plot_grid(pB,pC,ncol=1,rel_heights=c(1,1))
fig <- plot_grid(pA,right,ncol=2,rel_widths=c(1.32,1))

out <- "E:/222基金申请/2026年省自然/研究基础图6_TF筛选锁定FOXG1"
ggsave(paste0(out,".pdf"), fig, width=11.4, height=5.4, device=cairo_pdf)
ggsave(paste0(out,".png"), fig, width=11.4, height=5.4, dpi=300)
cat("SAVED:", out, ".pdf/.png\n")
