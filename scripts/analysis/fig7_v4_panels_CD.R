#!/usr/bin/env Rscript
# 图7 v4: C 改为通路级 z-score 热图(2 行 × 4 列), D 保留 violin
# 解决: 原 17 基因 z-score 不直观(snRNA 死亡基因表达噪声大), 改用 module score
A <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/anoikis_out"
OUT <- A
cell <- read.csv(file.path(A,"anoikis_mouse_per_cell.csv"))
cell$timepoint <- factor(cell$timepoint, levels=c("NoTBI","D1","D7","D28"))

navy<-"#1a3a5c"; mid<-"#3f6896"; steel<-"#7795b5"; accent<-"#b5563f"; teal<-"#2a6a5e"
gray1<-"#9aa7b3"; ink<-"#2a3a4a"

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

# ===== 计算 3x4 通路 z-score 矩阵 (加 anoikis 指数综合行) =====
mat <- matrix(0, 3, 4, dimnames=list(
  c("Death\n凋亡执行","Survival\nFAK/Akt 存活","Anoikis Index\n(Death − Survival)"),
  c("NoTBI","D1","D7","D28")))
death_mean   <- tapply(cell$death,    cell$timepoint, mean)
surv_mean    <- tapply(cell$survival, cell$timepoint, mean)
anoikis_mean <- tapply(cell$anoikis,  cell$timepoint, mean)
mat[1,] <- death_mean[c("NoTBI","D1","D7","D28")]
mat[2,] <- surv_mean[c("NoTBI","D1","D7","D28")]
mat[3,] <- anoikis_mean[c("NoTBI","D1","D7","D28")]

# 按行 z-score (跨时间点)
mat_z <- t(apply(mat, 1, scale))
colnames(mat_z) <- colnames(mat)

# 也计算 vs NoTBI 的 ΔScore (用于数值标注更直观)
fc_death   <- death_mean   - death_mean["NoTBI"]
fc_surv    <- surv_mean    - surv_mean["NoTBI"]
fc_anoikis <- anoikis_mean - anoikis_mean["NoTBI"]
fc_mat <- rbind(fc_death, fc_surv, fc_anoikis)

# 显著性 (vs NoTBI)
p_mat <- matrix(NA, 3, 4)
for (j in 1:4){
  tp <- colnames(mat)[j]
  if (tp == "NoTBI") { p_mat[,j] <- NA; next }
  base_d <- cell$death[cell$timepoint=="NoTBI"]
  v_d    <- cell$death[cell$timepoint==tp]
  p_mat[1,j] <- wilcox.test(v_d, base_d)$p.value
  base_s <- cell$survival[cell$timepoint=="NoTBI"]
  v_s    <- cell$survival[cell$timepoint==tp]
  p_mat[2,j] <- wilcox.test(v_s, base_s)$p.value
  base_a <- cell$anoikis[cell$timepoint=="NoTBI"]
  v_a    <- cell$anoikis[cell$timepoint==tp]
  p_mat[3,j] <- wilcox.test(v_a, base_a)$p.value
}

cat("=== 通路 mean (raw module score) ===\n"); print(round(mat, 4))
cat("\n=== 通路 z-score (跨时间点) ===\n"); print(round(mat_z, 3))
cat("\n=== 通路 ΔScore vs NoTBI ===\n"); print(round(fc_mat, 4))
cat("\n=== p values vs NoTBI ===\n"); print(p_mat)

# ===== 绘图 =====
draw <- function(){
  layout(matrix(c(1,2), 1, 2, byrow=TRUE), widths=c(1.25, 1))
  par(mar=c(5.5,9.5,4.6,5.0), mgp=c(2.8,0.7,0),
      family="Microsoft YaHei", cex.axis=0.98, cex.lab=1.05)

  # ===== Panel C: 2x4 通路 z-score heatmap =====
  nr <- nrow(mat_z); nc <- ncol(mat_z)
  zlim <- max(abs(mat_z))
  brk <- seq(-zlim, zlim, length.out=101)
  cols <- colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#F7F7F7","#F4A582","#D6604D","#B2182B"))(100)
  # 调高对比度: 用更陡渐变
  # image expects z[i,j] at (x[i], y[j]) — 我们要 x=时间点, y=通路 → 转置
  image(1:nc, 1:nr, t(mat_z)[, nr:1], col=cols, breaks=brk, axes=FALSE, xlab="", ylab="",
        useRaster=FALSE)
  rect(0.5, 0.5, nc+0.5, nr+0.5, border="black", lwd=2)
  # 分隔线 (三行之间)
  abline(h=1.5, col="black", lwd=1.6)
  abline(h=2.5, col="black", lwd=1.6)

  # X 轴 - 时间点
  axis(1, at=1:nc, labels=colnames(mat_z), las=1, tick=FALSE, cex.axis=1.15, font=2, line=-0.3)
  # Y 轴 - 通路 (image 显示是底到顶, 即 i=1 在底, i=3 在顶)
  # 我们想顶=凋亡执行, 中=存活, 底=Anoikis 指数 → t(mat_z)[, nr:1] 已反序
  # nr:1 = 3:1 → 顶 i=3 对应 mat 第 1 行(Death), 底 i=1 对应 mat 第 3 行(Anoikis Index)
  axis(2, at=1, labels="Anoikis Index\n(Death−Surv)",
       las=2, tick=FALSE, cex.axis=0.92, col.axis="purple4", font=2, line=-0.3)
  axis(2, at=2, labels="FAK/Akt\n存活",
       las=2, tick=FALSE, cex.axis=0.92, col.axis=navy,   font=2, line=-0.3)
  axis(2, at=3, labels="凋亡执行\n(Death)",
       las=2, tick=FALSE, cex.axis=0.92, col.axis=accent, font=2, line=-0.3)

  # 数值标注 (z-score 值 + 显著性星号)
  z_t <- t(mat_z)[, nr:1]   # 与 image 一致
  fc_t <- t(fc_mat)[, nr:1]
  p_t  <- t(p_mat)[, nr:1]
  for (j in 1:nc){
    for (i in 1:nr){
      z_val <- z_t[j, i]
      fc_val <- fc_t[j, i]
      p_val  <- p_t[j, i]
      txt_col <- if (abs(z_val) > zlim*0.55) "white" else ink
      # 在格子里写 ΔScore + 星号
      lab <- sprintf("%+.3f", fc_val)
      if (!is.na(p_val) && p_val < 0.05) lab <- paste0(lab, "\n", sig_label(p_val))
      text(j, i, lab, cex=0.85, col=txt_col, font=2)
    }
  }
  # 右侧 color bar
  par(xpd=NA)
  cb_x <- nc + 0.85; cb_w <- 0.18
  cb_y <- seq(0.5, nr+0.5, length.out=length(cols)+1)
  for (k in 1:length(cols))
    rect(cb_x, cb_y[k], cb_x+cb_w, cb_y[k+1], col=cols[k], border=NA)
  text(cb_x+cb_w/2, nr+0.85, "z-score", cex=0.78, font=2)
  text(cb_x+cb_w+0.05, c(0.5, (nr+1)/2, nr+0.5),
       c(sprintf("%+.1f",-zlim), "0", sprintf("%+.1f", zlim)),
       cex=0.78, adj=0)
  par(xpd=FALSE)

  # 字母 C 统一黑色加粗左上
  mtext("C", side=3, line=2.0, adj=-0.30, font=2, cex=2.2, col="black")
  # 格子值说明 (放底部 xlab 下方)
  mtext("(格内: ΔScore vs NoTBI + Wilcoxon p)", side=1, line=3.0, cex=0.78, col=ink)

  # ===== Panel D: 凋亡执行 per-cell violin =====
  par(mar=c(5.0,5.6,4.6,1.4))
  ylim <- range(cell$death, na.rm=TRUE)
  plot(NA, xlim=c(0.5,4.5), ylim=ylim*1.2, axes=FALSE,
       xlab="", ylab="凋亡执行通路评分 (per-cell)")
  abline(h=0, col=gray1, lty=2)
  axis(1, at=1:4, labels=levels(cell$timepoint), las=1, font=2, cex.axis=1.1)
  axis(2, cex.axis=0.95)
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
    text(i, y, sig_label(p), cex=1.6, col=accent, font=2)
  }
  mtext("D", side=3, line=2.0, adj=-0.20, font=2, cex=2.2, col="black")
}

png(file.path(OUT,"fig7_v4_panel_CD.png"), width=3000, height=1300, res=300,
    type="cairo", family="Microsoft YaHei")
draw(); dev.off()
cat("图7 v4 panels CD 已生成\n")
