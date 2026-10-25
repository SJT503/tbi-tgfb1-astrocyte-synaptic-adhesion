#!/usr/bin/env Rscript
# 统一研究基础图柱状图风格库
# 用法: source 此文件后调用 unified_bar() / unified_bar_dot()

# === 全局色板 ===
navy   <- "#1a3a5c"
mid    <- "#3f6896"
steel  <- "#7795b5"
accent <- "#b5563f"
accent_light <- "#d28876"
accent_dark  <- "#8d3a25"
teal   <- "#2a6a5e"
gray1  <- "#9aa7b3"
ink    <- "#2a3a4a"

# === 5 组配色 (Sham → TBI → TBI+TGFβ1 1/10/50ng) ===
PAL5 <- c(navy, steel, accent_light, accent, accent_dark)
# === 4 组配色 (Sham/CCI/CCI+SCR/CCI+KO) ===
PAL4 <- c(navy, accent, steel, teal)
# === 3 组配色 (Sham/Ipsi/Contra) ===
PAL3 <- c(navy, accent, steel)
# === 2 组 ===
PAL2 <- c(navy, accent)

sig_label <- function(p){
  if (is.na(p)) ""
  else if (p<1e-4) "****" else if (p<1e-3) "***"
  else if (p<0.01) "**" else if (p<0.05) "*" else "ns"
}

# 显著性括号绘制 (在 bar 上方)
draw_sig <- function(x1, x2, y, label, col=ink, lwd=1.4, cex=1.3, tick=0.8){
  segments(x1, y, x2, y, col=col, lwd=lwd)
  segments(x1, y, x1, y-tick, col=col, lwd=lwd)
  segments(x2, y, x2, y-tick, col=col, lwd=lwd)
  text((x1+x2)/2, y+tick*0.4, label, col=col, font=2, cex=cex)
}

# 统一柱状图 (bar+dot 或 bar only)
# means: 长度=n_group 的均值向量
# errs:  SE 向量
# labels: x 轴 group 名 (向量)
# cols: bar 填充色
# dots: 散点数据 list (每组 1 向量), NULL=不画散点
# sig_brackets: list( list(from=1,to=2,p=0.002) ... ) 显著性括号
unified_bar <- function(means, errs, labels, cols,
                        dots=NULL, sig_brackets=NULL,
                        ylab="", ylim=NULL, h0=TRUE,
                        bar_w=0.7, border_col=NA, border_lwd=0,
                        dot_pch=21, dot_cex=1.1, dot_lwd=1.2,
                        err_lwd=1.5, err_cap=0.04){
  n <- length(means)
  if (is.null(ylim)) {
    pad <- (max(means+errs, na.rm=TRUE) - min(c(means-errs, 0), na.rm=TRUE))*0.30
    ylim <- c(min(c(means-errs, 0), na.rm=TRUE) - pad*0.2,
              max(means+errs, na.rm=TRUE) + pad*1.0)
  }
  bp <- barplot(means, names.arg=NA, col=cols, border=border_col, ylim=ylim,
                ylab=ylab, axes=FALSE, cex.lab=1.05, width=bar_w, space=0.45)
  axis(2, cex.axis=0.95)
  if (h0) abline(h=0, col=gray1, lty=2, lwd=1)
  # error bars
  for (i in 1:n){
    arrows(bp[i], means[i]-errs[i], bp[i], means[i]+errs[i],
           angle=90, code=3, length=err_cap, col=ink, lwd=err_lwd)
  }
  # dots (jittered)
  if (!is.null(dots)){
    set.seed(42)
    for (i in 1:n){
      yv <- dots[[i]]
      if (length(yv)==0) next
      xv <- jitter(rep(bp[i], length(yv)), amount=0.18)
      points(xv, yv, pch=dot_pch, bg="white", col=ink,
             cex=dot_cex, lwd=dot_lwd)
    }
  }
  # x 轴 group label (倾斜)
  par(xpd=NA)
  text(bp, par("usr")[3] - (par("usr")[4]-par("usr")[3])*0.025,
       labels, srt=45, adj=c(1,1), cex=0.85, font=2, col=ink)
  par(xpd=FALSE)
  # 显著性括号
  if (!is.null(sig_brackets)){
    y_top <- max(means+errs, na.rm=TRUE)
    pad <- (par("usr")[4] - y_top) / max(1, length(sig_brackets) + 1)
    for (k in seq_along(sig_brackets)){
      sb <- sig_brackets[[k]]
      y_pos <- y_top + pad * k
      draw_sig(bp[sb$from], bp[sb$to], y_pos, sig_label(sb$p),
               tick=(par("usr")[4]-par("usr")[3])*0.012)
    }
  }
  invisible(bp)
}
