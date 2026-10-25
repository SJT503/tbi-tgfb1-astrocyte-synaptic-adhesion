#!/usr/bin/env Rscript
# Step 6b figure — neuronal module score change across TBI severity gradient
suppressMessages({library(ggplot2)})
OUT <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/step6_out"
d <- read.csv(file.path(OUT,"step6b_severity_gradient.csv"))

lab <- c(GSE290150="Moderate\n(GSE290150)", GSE271769="Severe\n(GSE271769)",
         GSE277487="Severe\n(GSE277487)", GSE209552="Human\n(GSE209552)")
d$dataset <- factor(lab[d$dataset], levels=lab)
d$module  <- factor(d$module, levels=c("Adhesion","Antiapoptotic","IntrinsicApoptosis"))
d$sig <- ifelse(d$pseudobulk_p<0.001,"***",
         ifelse(d$pseudobulk_p<0.01,"**",
         ifelse(d$pseudobulk_p<0.05,"*","ns")))
d$dir <- ifelse(d$delta>0,"up","down")

p <- ggplot(d, aes(dataset, delta, fill=dir)) +
  geom_col(width=0.66, colour="grey25", linewidth=0.3) +
  geom_hline(yintercept=0, linewidth=0.4, colour="grey40") +
  geom_text(aes(label=sig, vjust=ifelse(delta>=0,-0.3,1.2)), size=3, fontface="bold") +
  facet_wrap(~module, nrow=1) +
  scale_fill_manual(values=c(up="#E07B54", down="#4DAECD"), guide="none") +
  labs(title="Neuronal module score change (TBI - control) across severity gradient",
       subtitle="bars = mean module-score delta; asterisks = pseudobulk t-test (per-sample); ns = not significant",
       x=NULL, y="Module score delta (TBI - control)") +
  theme_bw(base_size=9) +
  theme(plot.title=element_text(face="bold",size=9.5),
        plot.subtitle=element_text(size=7.3,colour="grey35"),
        strip.background=element_rect(fill="grey92"),
        strip.text=element_text(face="bold",size=8),
        axis.text.x=element_text(size=7),
        panel.grid.minor=element_blank())
ggsave(file.path(OUT,"step6b_severity_gradient.pdf"), p, width=180, height=72, units="mm")
ggsave(file.path(OUT,"step6b_severity_gradient.png"), p, width=180, height=72, units="mm", dpi=300)
cat("saved step6b_severity_gradient figure\n")
