#!/usr/bin/env Rscript
# 图7 v2 (E/F 删减版): 仅 C heatmap + D 凋亡执行 violin
suppressMessages({})
A <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/anoikis_out"
OUT <- A

cell <- read.csv(file.path(A,"anoikis_mouse_per_cell.csv"))
gene <- read.csv(file.path(A,"anoikis_mouse_gene_expr.csv"))

cell$timepoint <- factor(cell$timepoint, levels=c("NoTBI","D1","D7","D28"))
gene$timepoint <- factor(gene$timepoint, levels=c("NoTBI","D1","D7","D28"))

navy<-"#1a3a5c"; mid<-"#3f6896"; steel<-"#7795b5"; accent<-"#b5563f"; teal<-"#2a6a5e"
gray1<-"#9aa7b3"; ink<-"#2a3a4a"

death_g <- intersect(c("Bcl2l11","Bmf","Bax","Bad","Casp9","Casp3","Cycs","Apaf1"), colnames(gene))
surv_g  <- intersect(c("Ptk2","Akt1","Akt2","Akt3","Bcl2","Mcl1","Bcl2l1","Birc5","Mapk1"), colnames(gene))
gene_order <- c(death_g, surv_g)

gm <- aggregate(gene[,gene_order], list(timepoint=gene$timepoint), mean)
rownames(gm) <- gm$timepoint; gm <- as.matrix(gm[,gene_order])
gm_z <- scale(gm); gm_z[is.na(gm_z)] <- 0

my_violin <- function(x_vec, group, at, col, width=0.4){
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
  par(mfrow=c(1,2), mar=c(4.8,4.8,3.6,1.4), mgp=c(2.8,0.7,0),
      family="Microsoft YaHei", cex.axis=0.95, cex.lab=1.05)

  # ===== Panel C: heatmap =====
  nr <- ncol(gm_z); nc <- nrow(gm_z)
  par(mar=c(5.5,5.4,3.6,2.5))
  zlim <- max(abs(gm_z))
  brk <- seq(-zlim, zlim, length.out=101)
  cols <- colorRampPalette(c("#2166AC","#67A9CF","#F7F7F7","#EF8A62","#B2182B"))(100)
  image(1:nc, 1:nr, gm_z, col=cols, breaks=brk, axes=FALSE, xlab="", ylab="", useRaster=FALSE)
  rect(0.5, 0.5, nc+0.5, nr+0.5, border="black", lwd=1.5)
  axis(1, at=1:nc, labels=rownames(gm_z), las=1, tick=FALSE, cex.axis=1.0, font=2)
  for (i in 1:nr){
    g <- gene_order[i]
    col_g <- if (g %in% death_g) accent else navy
    axis(2, at=i, labels=g, las=2, tick=FALSE, cex.axis=0.85, col.axis=col_g, font=2)
  }
  abline(h=length(death_g)+0.5, col="black", lwd=1.5, lty=2)
  text(nc+0.6, length(death_g)/2+0.5, "凋亡执行", srt=270, adj=0.5, cex=0.78, col=accent, font=2)
  text(nc+0.6, length(death_g)+length(surv_g)/2+0.5, "FAK/Akt 存活", srt=270, adj=0.5, cex=0.78, col=navy, font=2)
  mtext("C", side=3, line=1.6, adj=-0.20, font=2, cex=1.6)
  mtext("anoikis 核心基因表达 (z-score) 时序变化", side=3, line=0.4, cex=1.0, col=ink)

  # ===== Panel D: 凋亡执行 per-cell violin =====
  par(mar=c(4.8,5.4,3.6,1.4))
  ylim <- range(cell$death, na.rm=TRUE)
  plot(NA, xlim=c(0.5,4.5), ylim=ylim*1.2, axes=FALSE,
       xlab="", ylab="凋亡执行通路评分 (per-cell)")
  abline(h=0, col=gray1, lty=2)
  axis(1, at=1:4, labels=levels(cell$timepoint), las=1, font=2, cex.axis=1.0)
  axis(2)
  cols4 <- c(navy, accent, accent, accent)
  for (i in seq_along(levels(cell$timepoint))){
    tp <- levels(cell$timepoint)[i]
    x <- cell$death[cell$timepoint==tp]
    my_violin(x, tp, at=i, col=cols4[i])
  }
  base <- cell$death[cell$timepoint=="NoTBI"]
  for (i in 2:4){
    tp <- levels(cell$timepoint)[i]
    v <- cell$death[cell$timepoint==tp]
    p <- wilcox.test(v, base)$p.value
    y <- ylim[2]*1.05
    text(i, y, sig_label(p), cex=1.5, col=accent, font=2)
  }
  mtext("D", side=3, line=1.6, adj=-0.20, font=2, cex=1.6)
  mtext("凋亡执行通路 单细胞评分分布", side=3, line=0.4, cex=1.0, col=ink)
}

png(file.path(OUT,"fig7_v2_panel_CD.png"), width=2400, height=1200, res=300,
    type="cairo", family="Microsoft YaHei")
draw(); dev.off()
cat("图7 panels CD 已生成\n")
