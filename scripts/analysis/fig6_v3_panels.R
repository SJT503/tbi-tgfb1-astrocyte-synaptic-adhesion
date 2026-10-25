#!/usr/bin/env Rscript
# 图6 v3: B/C/D 升级为期刊感 heatmap+violin 风格
# A: TF screen scatter (单独保留, 复用现有)
# B: top TFs × 5 datasets heatmap (FOXG1 高亮)
# C: FOXG1 虚拟敲除对核心黏附基因的预测扰动 (cortex vs hippo dot+bar)
# D: 每细胞 FOXG1 表达 × 时间点 violin (替换原 per-sample bar)
suppressMessages({library(Seurat); library(Matrix)})
set.seed(1)
A <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis"
OUT <- file.path(A,"fig6_out")

# === 数据加载 ===
cross <- read.csv(file.path(A,"tf_screen_out","CROSS_3datasets_TFscreen.csv"))
g226  <- read.csv(file.path(A,"tf_screen_out","GSE226207_ms_TFscreen.csv"))
ginteg<- read.csv(file.path(A,"tf_screen_out","GSE101901_180862_ms_TFscreen.csv"))
g209  <- read.csv(file.path(A,"tf_screen_out","GSE209552_hs_TFscreen.csv"))
vko_ctx <- read.csv(file.path(A,"vko","foxg1_ko_shift.csv"))
vko_hip <- read.csv(file.path(A,"vko","foxg1_ko_hippo.csv"))
ps  <- read.csv(file.path(A,"foxg1_timecourse","foxg1_astro_per_sample.csv"))

# 5 个数据集 FOXG1 log2FC (与现有 v2 B 一致, 给 B 行参考)
foxg1_5ds <- data.frame(
  dataset = c("GSE180862\n海马","GSE180862\n皮层","GSE226207\n(纯度不合格)","GSE271769\n重度CCI","GSE209552\n人源"),
  lfc     = c(-0.28, -0.23, -0.24, -0.51, -0.76),
  rho     = c( 0.04,  0.08,  0.28,  0.20,  0.24)
)

# === 颜色 ===
navy<-"#1a3a5c"; mid<-"#3f6896"; steel<-"#7795b5"; accent<-"#b5563f"; teal<-"#2a6a5e"
gray1<-"#9aa7b3"; ink<-"#2a3a4a"

# ==== 工具 ====
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

# ==== 取每细胞 FOXG1 (用于 D) ====
rds <- file.path(A,"foxg1_timecourse","GSE271769_astro.rds")
cat("加载 GSE271769 astro RDS...\n")
seu <- readRDS(rds)
DefaultAssay(seu) <- "RNA"
# 纯度门禁 (此 RDS 已是 astro subset)
expr_check <- GetAssayData(seu, layer="data")
markers_keep <- c("Aqp4","Gfap","Snap25","Syt1")
present <- intersect(markers_keep, rownames(expr_check))
if (all(c("Aqp4","Snap25") %in% rownames(expr_check))) {
  aqp <- mean(expr_check["Aqp4",]); snap <- mean(expr_check["Snap25",])
  cat(sprintf("纯度 Aqp4/Snap25 = %.2f\n", aqp/max(snap,1e-6)))
  if (aqp/max(snap,1e-6) < 3) warning("Astro 纯度不足, 慎用")
}
# 提取 Foxg1
foxg1_vec <- as.numeric(GetAssayData(seu, layer="data")["Foxg1",])
tp_vec    <- as.character(seu$timepoint)
percell <- data.frame(foxg1=foxg1_vec, timepoint=tp_vec)
percell$timepoint <- factor(percell$timepoint, levels=c("NoTBI","D1","D7","D28"))
cat("每细胞 FOXG1 n=", nrow(percell), "\n")

# ==== 取 top down-regulated TFs (B 用) ====
# 按 mean_lfc 升序 + n_lostact >= 2 (黏附基因失活数)
top_tf <- cross[order(cross$mean_lfc), ]
top_tf <- top_tf[top_tf$n_lostact >= 2 & !is.na(top_tf$mean_lfc), ]
# 强制 FOXG1 入选
if (!"FOXG1" %in% head(top_tf$TF, 12)) {
  top_tf <- rbind(top_tf[top_tf$TF=="FOXG1",], top_tf[top_tf$TF!="FOXG1",])
}
top12 <- head(top_tf$TF, 12)
cat("Top 12 TF:", top12, "\n")

# 构建 TF × dataset log2FC 矩阵 (用 3 个 cross 列 + 5ds 数据中 FOXG1 单独标注)
mat_b <- as.matrix(cross[match(top12, cross$TF), c("lfc_integ","lfc_226207","lfc_209552")])
colnames(mat_b) <- c("Integ\n(GSE101901+180862)","GSE226207\n(纯度低)","GSE209552\n(人源)")
rownames(mat_b) <- top12
mat_b[is.na(mat_b)] <- 0

# ==== 取核心黏附基因 KO shift (C 用) ====
adh_genes <- c("Nrxn1","Nrxn2","Nrxn3","Cadm1","Cadm2","Cntn1","Cntn2","Cdh2","Nlgn1","Nlgn2")
ctx_sub <- vko_ctx[vko_ctx$gene %in% adh_genes, ]
hip_sub <- vko_hip[vko_hip$gene %in% adh_genes, ]
# 合并: 取交集顺序
adh_use <- intersect(adh_genes, ctx_sub$gene)
adh_use <- intersect(adh_use, hip_sub$gene)
cat("黏附基因 in vko:", adh_use, "\n")
ctx_shift <- ctx_sub$ko_shift[match(adh_use, ctx_sub$gene)]
hip_shift <- hip_sub$ko_shift[match(adh_use, hip_sub$gene)]

# ==== 总体绘图 ====
draw <- function(){
  layout(matrix(c(1,2,3,4), 2, 2, byrow=TRUE), heights=c(1,1), widths=c(1,1))
  par(mar=c(4.8,5.0,3.6,1.4), mgp=c(2.8,0.7,0),
      family="Microsoft YaHei", cex.axis=0.92, cex.lab=1.0)

  # ===== Panel A: TF screen scatter (重绘以风格统一) =====
  par(mar=c(4.8,5.0,3.6,1.4))
  x <- ginteg$lfc; y <- ginteg$rho_sa
  if (is.null(y)) y <- ginteg$rho
  if (is.null(y)) {
    # 试探备选列名
    yc <- grep("rho|cor", colnames(ginteg), value=TRUE, ignore.case=TRUE)
    cat("rho 列候选:", yc, "\n")
    y <- ginteg[[ yc[1] ]]
  }
  # 灰底散点 (其他 TF 弱化)
  plot(x, y, pch=21, bg=alpha_col("#9aa7b3",0.30), col=alpha_col("#9aa7b3",0.50), cex=0.75,
       xlab="损伤效应  log2FC  (TBI vs Sham)",
       ylab="与突触黏附程序的 Spearman r",
       xlim=c(min(x,na.rm=T)-0.1, max(x,na.rm=T)+0.15),
       ylim=c(min(y,na.rm=T)-0.02, max(y,na.rm=T)+0.06))
  abline(v=0, col=gray1, lty=2); abline(h=0, col=gray1, lty=2)
  # 象限淡背景填充 (左上=目标象限 浅色)
  rect(par("usr")[1], 0, 0, par("usr")[4],
       col=alpha_col(accent, 0.06), border=NA)
  # 象限提示
  text(par("usr")[1]+0.05, par("usr")[4]-0.005,
       "损伤下调 × 与黏附正相关\n→ 丢失的黏附正调控者",
       cex=0.80, col=ink, adj=c(0,1), font=2)
  # FOXG1 突出: 同心圆 (光晕) + 大点 + 引线 + 大字标签
  # 兼容大小写: 小鼠数据集 TF 列为 "Foxg1", 用 TF_upper 或 ignore.case 匹配
  if ("TF_upper" %in% colnames(ginteg)) {
    i <- which(ginteg$TF_upper == "FOXG1")
  } else {
    i <- which(toupper(ginteg$TF) == "FOXG1")
  }
  if (length(i)) {
    fx <- x[i]; fy <- y[i]
    # 光晕圈 (三层渐变)
    points(fx, fy, pch=21, bg=NA, col=alpha_col(accent,0.20), cex=4.8, lwd=12)
    points(fx, fy, pch=21, bg=NA, col=alpha_col(accent,0.40), cex=3.5, lwd=8)
    points(fx, fy, pch=21, bg=alpha_col(accent,0.95), col=accent, cex=2.2, lwd=2)
    # 引线 + 标签
    lx <- fx + 0.45; ly <- fy + 0.10
    segments(fx+0.06, fy+0.005, lx-0.02, ly-0.005, col=accent, lwd=1.6)
    rect(lx-0.02, ly-0.011, lx+0.55, ly+0.018, col="white", border=accent, lwd=1.6)
    text(lx+0.02, ly+0.003, "FOXG1", col=accent, font=2, cex=1.4, adj=0)
  }
  mtext("A", side=3, line=1.6, adj=-0.20, font=2, cex=2.2, col="black")

  # ===== Panel B: top 12 TFs × 3 datasets heatmap =====
  par(mar=c(5.8,5.0,3.6,3.0))
  nr <- nrow(mat_b); nc <- ncol(mat_b)
  zlim <- max(abs(mat_b))
  brk <- seq(-zlim, zlim, length.out=101)
  cols <- colorRampPalette(c("#2166AC","#67A9CF","#F7F7F7","#EF8A62","#B2182B"))(100)
  image(1:nc, 1:nr, t(mat_b)[, nr:1], col=cols, breaks=brk, axes=FALSE, xlab="", ylab="")
  rect(0.5, 0.5, nc+0.5, nr+0.5, border="black", lwd=1.5)
  axis(1, at=1:nc, labels=colnames(mat_b), las=1, tick=FALSE, cex.axis=0.78, font=2, line=-0.2)
  for (i in 1:nr){
    g <- top12[nr - i + 1]
    col_g <- if (g == "FOXG1") accent else ink
    fnt   <- if (g == "FOXG1") 2 else 1
    axis(2, at=i, labels=g, las=2, tick=FALSE, cex.axis=0.82, col.axis=col_g, font=fnt)
    # 数值标注
    for (j in 1:nc){
      val <- t(mat_b)[, nr:1][j, i]
      txt_col <- if (abs(val) > zlim*0.5) "white" else ink
      text(j, i, sprintf("%.2f", val), cex=0.62, col=txt_col, font=if(g=="FOXG1") 2 else 1)
    }
  }
  # FOXG1 行框高亮
  i_foxg1 <- which(top12 == "FOXG1")
  if (length(i_foxg1)){
    i_plot <- nr - i_foxg1 + 1
    rect(0.5, i_plot-0.5, nc+0.5, i_plot+0.5, border=accent, lwd=2.4, lty=1)
  }
  # color bar
  par(xpd=NA)
  cb_x <- nc + 0.85; cb_w <- 0.18
  cb_y <- seq(0.5, nr+0.5, length.out=length(cols)+1)
  for (k in 1:length(cols))
    rect(cb_x, cb_y[k], cb_x+cb_w, cb_y[k+1], col=cols[k], border=NA)
  text(cb_x+cb_w/2, nr+0.9, "log2FC", cex=0.7, font=2)
  text(cb_x+cb_w+0.05, c(0.5, (nr+1)/2, nr+0.5),
       c(sprintf("%.1f",-zlim), "0", sprintf("%.1f",zlim)),
       cex=0.65, adj=0)
  par(xpd=FALSE)
  mtext("B", side=3, line=1.6, adj=-0.20, font=2, cex=2.2, col="black")

  # ===== Panel C: FOXG1-KO 黏附基因扰动 (cortex vs hippo) =====
  # 颜色: 皮层=深海蓝 navy, 海马=橙红 accent (对比度大)
  par(mar=c(7.8,5.4,3.6,1.6))   # 加大底部 mar 容纳图例
  nn <- length(adh_use)
  ctx_col <- navy
  hip_col <- accent
  x_lim <- range(c(ctx_shift, hip_shift), na.rm=TRUE)
  pad <- (max(x_lim)-min(x_lim))*0.12
  x_lim <- c(min(x_lim)-pad, max(0.005, max(x_lim))+pad)
  plot(NA, xlim=x_lim, ylim=c(0.5, nn+0.5), axes=FALSE,
       xlab="预测 KO 后基因表达扰动 (mean_shift)", ylab="")
  abline(v=0, col=gray1, lty=2)
  # 平行 lollipop / dumbbell
  for (i in 1:nn){
    y0 <- i
    # cortex (navy, 上)
    segments(0, y0+0.20, ctx_shift[i], y0+0.20, col=ctx_col, lwd=3)
    points(ctx_shift[i], y0+0.20, pch=21, bg=ctx_col, col=ctx_col, cex=1.5)
    # hippo (accent, 下)
    segments(0, y0-0.20, hip_shift[i], y0-0.20, col=hip_col, lwd=3)
    points(hip_shift[i], y0-0.20, pch=23, bg=hip_col, col=hip_col, cex=1.5)
  }
  axis(2, at=1:nn, labels=adh_use, las=2, tick=FALSE, cex.axis=0.88, font=2)
  axis(1, cex.axis=0.92)
  # 图例: 移到底部 mar 区 (绘图区外, 不遮挡)
  par(xpd=NA)
  leg_y <- par("usr")[3] - (par("usr")[4]-par("usr")[3]) * 0.18
  leg_x <- (par("usr")[1] + par("usr")[2]) / 2
  legend(leg_x, leg_y, xjust=0.5, yjust=0.5,
         legend=c("皮层 GSE271769", "海马 整合 (GSE101901+180862)"),
         pch=c(21, 23), pt.bg=c(ctx_col, hip_col), col=c(ctx_col, hip_col),
         bty="n", cex=0.85, pt.cex=1.4, horiz=TRUE, text.width=NA)
  par(xpd=FALSE)
  mtext("C", side=3, line=1.6, adj=-0.20, font=2, cex=2.2, col="black")

  # ===== Panel D: per-cell FOXG1 violin × 时间点 =====
  par(mar=c(4.8,5.0,3.6,1.4))
  ylim <- range(percell$foxg1, na.rm=TRUE)
  plot(NA, xlim=c(0.5,4.5), ylim=c(ylim[1], ylim[2]*1.15), axes=FALSE,
       xlab="", ylab="FOXG1 表达 (log-normalized, per-cell)")
  abline(h=median(percell$foxg1[percell$timepoint=="NoTBI"], na.rm=TRUE),
         col=gray1, lty=2)
  axis(1, at=1:4, labels=levels(percell$timepoint), las=1, font=2)
  axis(2)
  cols4 <- c(navy, accent, mid, accent)
  for (i in seq_along(levels(percell$timepoint))){
    tp <- levels(percell$timepoint)[i]
    xv <- percell$foxg1[percell$timepoint==tp]
    my_violin(xv, at=i, col=cols4[i])
  }
  base <- percell$foxg1[percell$timepoint=="NoTBI"]
  for (i in 2:4){
    tp <- levels(percell$timepoint)[i]
    v <- percell$foxg1[percell$timepoint==tp]
    p <- wilcox.test(v, base)$p.value
    y <- ylim[2]*1.05
    text(i, y, sig_label(p), cex=1.4, col=cols4[i], font=2)
  }
  mtext("D", side=3, line=1.6, adj=-0.20, font=2, cex=2.2, col="black")
}

# 辅助 alpha (R 4 内置: adjustcolor)
alpha_col <- function(c, a) adjustcolor(c, alpha.f=a)

png(file.path(OUT,"fig6_v3_panel_ABCD.png"), width=2800, height=2400, res=300,
    type="cairo", family="Microsoft YaHei")
draw(); dev.off()
cat("图6 v3 已生成\n")
