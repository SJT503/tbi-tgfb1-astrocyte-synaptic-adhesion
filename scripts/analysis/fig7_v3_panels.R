#!/usr/bin/env Rscript
# 图7 v3: 6 面板 (按用户原指令: A/B TUNEL + C/D 小鼠 + E/F 人源)
# C: 17 anoikis 基因 heatmap × 4 时间点 (小鼠 GSE271769)
# D: 凋亡执行通路 per-cell violin (小鼠 4 时间点)
# E: 凋亡执行通路 per-cell violin (人源 GSE209552 CTL vs TBI)
# F: 凋亡执行通路 per-sample bar+dots (人源 CTL n=5 vs TBI n=12)
suppressMessages({})
A <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/anoikis_out"
OUT <- A

cell <- read.csv(file.path(A,"anoikis_mouse_per_cell.csv"))
gene <- read.csv(file.path(A,"anoikis_mouse_gene_expr.csv"))
hum_c <- read.csv(file.path(A,"anoikis_human_per_cell.csv"))
hum_s <- read.csv(file.path(A,"anoikis_human_per_sample.csv"))

cell$timepoint <- factor(cell$timepoint, levels=c("NoTBI","D1","D7","D28"))
gene$timepoint <- factor(gene$timepoint, levels=c("NoTBI","D1","D7","D28"))
hum_c$group <- factor(hum_c$group, levels=c("CTL","TBI"))
hum_s$group <- factor(hum_s$group, levels=c("CTL","TBI"))

navy<-"#1a3a5c"; mid<-"#3f6896"; steel<-"#7795b5"; accent<-"#b5563f"; teal<-"#2a6a5e"
gray1<-"#9aa7b3"; ink<-"#2a3a4a"

death_g <- intersect(c("Bcl2l11","Bmf","Bax","Bad","Casp9","Casp3","Cycs","Apaf1"), colnames(gene))
surv_g  <- intersect(c("Ptk2","Akt1","Akt2","Akt3","Bcl2","Mcl1","Bcl2l1","Birc5","Mapk1"), colnames(gene))
gene_order <- c(death_g, surv_g)

gm <- aggregate(gene[,gene_order], list(timepoint=gene$timepoint), mean)
rownames(gm) <- gm$timepoint; gm <- as.matrix(gm[,gene_order])
gm_z <- scale(gm); gm_z[is.na(gm_z)] <- 0

my_violin <- function(x_vec, at, col, width=0.4){
  d <- density(x_vec, na.rm=TRUE, n=128)
  d$y <- d$y / max(d$y) * width
  polygon(c(at-d$y, rev(at+d$y)), c(d$x, rev(d$x)), col=col, border=col)
  med <- median(x_vec, na.rm=TRUE)
  segments(at-width*0.6, med, at+width*0.6, med, col="white", lwd=2.5)
  q <- quantile(x_vec, c(0.25, 0.75), na.rm=TRUE)
  segments(at, q[1], at, q[2], col="white", lwd=4)
}
sig_label <- function(p){
  if (p<1e-4) "****" else if (p<1e-3) "***" else if (p<0.01) "**" else if (p<0.05) "*" else "ns"
}

draw <- function(){
  par(mfrow=c(2,2), mar=c(4.8,4.8,3.6,1.4), mgp=c(2.8,0.7,0),
      family="Microsoft YaHei", cex.axis=0.92, cex.lab=1.0)

  # ==== C: 小鼠 heatmap ====
  nr <- ncol(gm_z); nc <- nrow(gm_z)
  par(mar=c(5.5,5.0,3.6,2.5))
  zlim <- max(abs(gm_z)); brk <- seq(-zlim, zlim, length.out=101)
  cols <- colorRampPalette(c("#2166AC","#67A9CF","#F7F7F7","#EF8A62","#B2182B"))(100)
  image(1:nc, 1:nr, gm_z, col=cols, breaks=brk, axes=FALSE, xlab="", ylab="", useRaster=FALSE)
  rect(0.5, 0.5, nc+0.5, nr+0.5, border="black", lwd=1.5)
  axis(1, at=1:nc, labels=rownames(gm_z), las=1, tick=FALSE, cex.axis=0.95, font=2)
  for (i in 1:nr){
    g <- gene_order[i]
    col_g <- if (g %in% death_g) accent else navy
    axis(2, at=i, labels=g, las=2, tick=FALSE, cex.axis=0.78, col.axis=col_g, font=2)
  }
  abline(h=length(death_g)+0.5, col="black", lwd=1.5, lty=2)
  text(nc+0.6, length(death_g)/2+0.5, "凋亡执行", srt=270, adj=0.5, cex=0.7, col=accent, font=2)
  text(nc+0.6, length(death_g)+length(surv_g)/2+0.5, "FAK/Akt 存活", srt=270, adj=0.5, cex=0.7, col=navy, font=2)
  mtext("C", side=3, line=1.6, adj=-0.20, font=2, cex=1.5)
  mtext("小鼠 GSE271769  anoikis 核心基因 z-score 时序", side=3, line=0.4, cex=0.88, col=ink)

  # ==== D: 小鼠 per-cell death violin ====
  par(mar=c(4.8,5.0,3.6,1.4))
  ylim <- range(cell$death, na.rm=TRUE)
  plot(NA, xlim=c(0.5,4.5), ylim=ylim*1.2, axes=FALSE,
       xlab="", ylab="凋亡执行通路评分 (per-cell)")
  abline(h=0, col=gray1, lty=2)
  axis(1, at=1:4, labels=levels(cell$timepoint), las=1, font=2)
  axis(2)
  cols4 <- c(navy, accent, accent, accent)
  for (i in seq_along(levels(cell$timepoint))){
    tp <- levels(cell$timepoint)[i]
    x <- cell$death[cell$timepoint==tp]
    my_violin(x, at=i, col=cols4[i])
  }
  base <- cell$death[cell$timepoint=="NoTBI"]
  for (i in 2:4){
    tp <- levels(cell$timepoint)[i]
    v <- cell$death[cell$timepoint==tp]
    p <- wilcox.test(v, base)$p.value
    y <- ylim[2]*1.05
    text(i, y, sig_label(p), cex=1.4, col=accent, font=2)
  }
  mtext("D", side=3, line=1.6, adj=-0.20, font=2, cex=1.5)
  mtext("小鼠重度 CCI 凋亡执行通路单细胞分布", side=3, line=0.4, cex=0.88, col=ink)

  # ==== E: 人源 per-cell death violin (CTL vs TBI) ====
  par(mar=c(4.8,5.0,3.6,1.4))
  ylim <- range(hum_c$death, na.rm=TRUE)
  plot(NA, xlim=c(0.5,2.5), ylim=ylim*1.2, axes=FALSE,
       xlab="", ylab="凋亡执行通路评分 (per-cell)")
  abline(h=0, col=gray1, lty=2)
  axis(1, at=1:2, labels=c("CTL","TBI"), las=1, font=2)
  axis(2)
  cols2 <- c(navy, accent)
  for (i in seq_along(levels(hum_c$group))){
    gp <- levels(hum_c$group)[i]
    x <- hum_c$death[hum_c$group==gp]
    my_violin(x, at=i, col=cols2[i])
  }
  base <- hum_c$death[hum_c$group=="CTL"]
  v    <- hum_c$death[hum_c$group=="TBI"]
  p_e  <- wilcox.test(v, base)$p.value
  p_str <- if (p_e < 1e-4) sprintf("p<1e-4 %s (per-cell n 极大)", sig_label(p_e)) else sprintf("p=%.3g %s", p_e, sig_label(p_e))
  text(1.5, ylim[2]*1.05, p_str, cex=0.85, col=ink, font=2)
  n_ctl <- sum(hum_c$group=="CTL"); n_tbi <- sum(hum_c$group=="TBI")
  text(1, ylim[1]*1.15, sprintf("n=%d cells", n_ctl), cex=0.75, col=ink)
  text(2, ylim[1]*1.15, sprintf("n=%d cells", n_tbi), cex=0.75, col=ink)
  mtext("E", side=3, line=1.6, adj=-0.20, font=2, cex=1.5)
  mtext("人源 GSE209552 神经元单细胞分布 (慢性期幸存)", side=3, line=0.4, cex=0.88, col=ink)

  # ==== F: 人源 per-sample bar (CTL n=5 vs TBI n=12) ====
  par(mar=c(4.8,5.0,3.6,1.4))
  m  <- tapply(hum_s$death, hum_s$group, mean)
  se <- tapply(hum_s$death, hum_s$group, function(x) sd(x)/sqrt(length(x)))
  ymax <- max(m+se); ymin <- min(c(0, m-se))
  pad <- (ymax-ymin)*0.30
  ylim <- c(ymin-pad*0.3, ymax+pad)
  bp <- barplot(m, col=c(navy,accent), border=NA, ylim=ylim,
                ylab="凋亡执行通路评分 (per-sample)",
                names.arg=names(m), cex.names=0.95, font.axis=2)
  abline(h=0, col=gray1, lty=2)
  arrows(bp, m-se, bp, m+se, angle=90, code=3, length=0.05, col=ink, lwd=1.4)
  set.seed(1)
  for (i in seq_along(levels(hum_s$group))){
    gp <- levels(hum_s$group)[i]
    y <- hum_s$death[hum_s$group==gp]
    points(jitter(rep(bp[i],length(y)),amount=0.12), y, pch=21, bg="white", col=ink, cex=1.1, lwd=1.2)
  }
  base <- hum_s$death[hum_s$group=="CTL"]
  v    <- hum_s$death[hum_s$group=="TBI"]
  p_f  <- wilcox.test(v, base)$p.value
  text(mean(bp), ymax+pad*0.6, sprintf("p=%.2f %s", p_f, sig_label(p_f)), cex=0.95, col=ink, font=2)
  text(bp[1], ymin-pad*0.15, "n=5", cex=0.78, col=ink)
  text(bp[2], ymin-pad*0.15, "n=12", cex=0.78, col=ink)
  mtext("F", side=3, line=1.6, adj=-0.20, font=2, cex=1.5)
  mtext("人源 GSE209552 生物学重复 (慢性期幸存)", side=3, line=0.4, cex=0.88, col=ink)
}

png(file.path(OUT,"fig7_v3_panel_CDEF.png"), width=2400, height=2400, res=300,
    type="cairo", family="Microsoft YaHei")
draw(); dev.off()
cat("图7 v3 panels CDEF 已生成 (C/D 小鼠 + E/F 人源)\n")
