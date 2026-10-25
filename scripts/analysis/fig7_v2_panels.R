#!/usr/bin/env Rscript
# 图7 v2: 小鼠 GSE271769 深化版统计面板 (C/D/E/F)
# C: 17 anoikis 基因 heatmap × 4 时间点
# D: 凋亡执行通路 per-cell violin
# E: 存活通路 per-cell violin
# F: 失巢凋亡指数 per-sample bar
suppressMessages({})
A <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/anoikis_out"
OUT <- A

# 加载
cell <- read.csv(file.path(A,"anoikis_mouse_per_cell.csv"))
samp <- read.csv(file.path(A,"anoikis_per_sample.csv"))
gene <- read.csv(file.path(A,"anoikis_mouse_gene_expr.csv"))

cell$timepoint <- factor(cell$timepoint, levels=c("NoTBI","D1","D7","D28"))
samp$timepoint <- factor(samp$timepoint, levels=c("NoTBI","D1","D7","D28"))
gene$timepoint <- factor(gene$timepoint, levels=c("NoTBI","D1","D7","D28"))

navy<-"#1a3a5c"; mid<-"#3f6896"; steel<-"#7795b5"; accent<-"#b5563f"; teal<-"#2a6a5e"
gray1<-"#9aa7b3"; ink<-"#2a3a4a"

# 基因分组(死亡 / 存活)
death_g <- intersect(c("Bcl2l11","Bmf","Bax","Bad","Casp9","Casp3","Cycs","Apaf1"), colnames(gene))
surv_g  <- intersect(c("Ptk2","Akt1","Akt2","Akt3","Bcl2","Mcl1","Bcl2l1","Birc5","Mapk1"), colnames(gene))
gene_order <- c(death_g, surv_g)

# 每基因每时间点 mean
gm <- aggregate(gene[,gene_order], list(timepoint=gene$timepoint), mean)
rownames(gm) <- gm$timepoint; gm <- as.matrix(gm[,gene_order])
# z-score 每基因 (列方向)
gm_z <- scale(gm)  # 默认列z-score
gm_z[is.na(gm_z)] <- 0

# 自定义 violin (无依赖)
my_violin <- function(x_vec, group, at, col, width=0.4){
  d <- density(x_vec, na.rm=TRUE, n=128)
  d$y <- d$y / max(d$y) * width
  polygon(c(at-d$y, rev(at+d$y)), c(d$x, rev(d$x)), col=col, border=col)
  # 中位数横线
  med <- median(x_vec, na.rm=TRUE)
  segments(at-width*0.6, med, at+width*0.6, med, col="white", lwd=2.5)
  # IQR
  q <- quantile(x_vec, c(0.25, 0.75), na.rm=TRUE)
  segments(at, q[1], at, q[2], col="white", lwd=4)
}

# 显著性 (Mann-Whitney 与 NoTBI 比)
sig_label <- function(p){
  if (p<1e-4) "****" else if (p<1e-3) "***" else if (p<0.01) "**" else if (p<0.05) "*" else "ns"
}

draw <- function(){
  par(mfrow=c(2,2), mar=c(4.8,4.8,3.6,1.4), mgp=c(2.8,0.7,0),
      family="Microsoft YaHei", cex.axis=0.92, cex.lab=1.0)

  # ========== Panel C: heatmap ==========
  nr <- ncol(gm_z); nc <- nrow(gm_z)
  par(mar=c(5.5,5.0,3.6,2.5))
  # 颜色: blue-white-red
  zlim <- max(abs(gm_z))
  brk <- seq(-zlim, zlim, length.out=101)
  cols <- colorRampPalette(c("#2166AC","#67A9CF","#F7F7F7","#EF8A62","#B2182B"))(100)
  image(1:nc, 1:nr, gm_z, col=cols, breaks=brk, axes=FALSE, xlab="", ylab="",
        useRaster=FALSE)
  # 边框
  rect(0.5, 0.5, nc+0.5, nr+0.5, border="black", lwd=1.5)
  # 列轴 (时间点)
  axis(1, at=1:nc, labels=rownames(gm_z), las=1, tick=FALSE, cex.axis=0.95, font=2)
  # 行轴 (基因) — 死亡基因红色, 存活基因蓝色
  for (i in 1:nr){
    g <- gene_order[i]
    col_g <- if (g %in% death_g) accent else navy
    axis(2, at=i, labels=g, las=2, tick=FALSE, cex.axis=0.78, col.axis=col_g, font=2)
  }
  # 右侧分组标识 (death vs survival)
  abline(h=length(death_g)+0.5, col="black", lwd=1.5, lty=2)
  text(nc+0.6, length(death_g)/2+0.5, "凋亡执行", srt=270, adj=0.5, cex=0.7, col=accent, font=2)
  text(nc+0.6, length(death_g)+length(surv_g)/2+0.5, "FAK/Akt 存活", srt=270, adj=0.5, cex=0.7, col=navy, font=2)
  mtext("C", side=3, line=1.6, adj=-0.20, font=2, cex=1.5)
  mtext("anoikis 核心基因表达 (z-score) 时序变化", side=3, line=0.4, cex=0.93, col=ink)

  # ========== Panel D: 凋亡执行通路 per-cell violin ==========
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
    my_violin(x, tp, at=i, col=cols4[i])
  }
  # 显著性
  base <- cell$death[cell$timepoint=="NoTBI"]
  for (i in 2:4){
    tp <- levels(cell$timepoint)[i]
    v <- cell$death[cell$timepoint==tp]
    p <- wilcox.test(v, base)$p.value
    y <- ylim[2]*1.05
    text(i, y, sig_label(p), cex=1.4, col=accent, font=2)
  }
  mtext("D", side=3, line=1.6, adj=-0.20, font=2, cex=1.5)
  mtext("凋亡执行通路 单细胞评分分布", side=3, line=0.4, cex=0.93, col=ink)

  # ========== Panel E: 存活通路 per-cell violin ==========
  par(mar=c(4.8,5.0,3.6,1.4))
  ylim <- range(cell$survival, na.rm=TRUE)
  plot(NA, xlim=c(0.5,4.5), ylim=ylim*1.2, axes=FALSE,
       xlab="", ylab="FAK/Akt 存活通路评分 (per-cell)")
  abline(h=0, col=gray1, lty=2)
  axis(1, at=1:4, labels=levels(cell$timepoint), las=1, font=2)
  axis(2)
  cols4e <- c(navy, mid, mid, mid)
  for (i in seq_along(levels(cell$timepoint))){
    tp <- levels(cell$timepoint)[i]
    x <- cell$survival[cell$timepoint==tp]
    my_violin(x, tp, at=i, col=cols4e[i])
  }
  base <- cell$survival[cell$timepoint=="NoTBI"]
  for (i in 2:4){
    tp <- levels(cell$timepoint)[i]
    v <- cell$survival[cell$timepoint==tp]
    p <- wilcox.test(v, base)$p.value
    y <- ylim[2]*1.05
    text(i, y, sig_label(p), cex=1.4, col=mid, font=2)
  }
  mtext("E", side=3, line=1.6, adj=-0.20, font=2, cex=1.5)
  mtext("FAK/Akt 存活通路 单细胞评分分布", side=3, line=0.4, cex=0.93, col=ink)

  # ========== Panel F: 失巢凋亡指数 per-sample bar + dots ==========
  par(mar=c(4.8,5.0,3.6,1.4))
  # 直接用死亡通路评分 (death module) — 跟 anoikis 指数同向但更直观:
  # 数值上升 = 凋亡执行通路升高 = anoikis 倾向升高
  m  <- tapply(samp$death, samp$timepoint, mean)
  se <- tapply(samp$death, samp$timepoint, function(x) sd(x)/sqrt(length(x)))
  ymax <- max(m+se); ymin <- min(m-se)
  pad <- (ymax-ymin)*0.30
  ylim <- c(ymin-pad*0.3, ymax+pad)
  bp <- barplot(m, col=c(navy,steel,accent,accent), border=NA, ylim=ylim,
                ylab="凋亡执行通路评分  (per-sample mean)",
                main="", names.arg=names(m), cex.names=0.95, font.axis=2)
  abline(h=0, col=gray1, lty=2)
  arrows(bp, m-se, bp, m+se, angle=90, code=3, length=0.05, col=ink, lwd=1.4)
  set.seed(1)
  for (i in seq_along(levels(samp$timepoint))){
    tp <- levels(samp$timepoint)[i]; y <- samp$death[samp$timepoint==tp]
    points(jitter(rep(bp[i],length(y)),amount=0.10), y, pch=21, bg="white", col=ink, cex=1.1, lwd=1.2)
  }
  base <- samp$death[samp$timepoint=="NoTBI"]
  for (i in 2:4){
    tp <- levels(samp$timepoint)[i]
    v <- samp$death[samp$timepoint==tp]
    p <- wilcox.test(v, base)$p.value
    if (p < 0.05){
      text(bp[i], m[i]+se[i]+pad*0.25, sig_label(p), cex=1.4, col=accent, font=2)
    }
  }
  mtext("F", side=3, line=1.6, adj=-0.20, font=2, cex=1.5)
  mtext("凋亡执行通路 生物学重复（n=4/组）", side=3, line=0.4, cex=0.93, col=ink)
}

png(file.path(OUT,"fig7_v2_panel_CDEF.png"), width=2400, height=2400, res=300,
    type="cairo", family="Microsoft YaHei")
draw(); dev.off()
cat("图7 v2 panels CDEF 已生成\n")
