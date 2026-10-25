#!/usr/bin/env Rscript
# 研究基础新增图 v3 (3 面板, 多度量综合证据):
# A: NGL 通路 3 组 (强支持)
# B: NRXN 多度量 3 组 (LR平均强度 + Astro_C1 输出)
# C: 亚群 × 通路 log2FC 热图 (异质性透明展示)
suppressMessages({library(CellChat)})
source("D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/style_lib.R")

D <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/wet_lab_data/results/cellchat_nrxn"
OUT <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/fig_selftest_v3.png"

cc <- list(
  Sham = readRDS(file.path(D,"cellchat_nrxn_Sham.rds")),
  CCI  = readRDS(file.path(D,"cellchat_nrxn_CCI_Scr.rds")),
  KO   = readRDS(file.path(D,"cellchat_nrxn_CCI_TGFbetaR1_KO.rds"))
)
ASTROS <- c("Astro_C0","Astro_C1","Astro_C2")
NEURONS <- c("Excitatory_Neuron","Inhibitory_Neuron")

# ===== 数据 1 (Panel A): NGL Astro→Exc_Neuron 总强度 (sum) =====
ngl_sum <- numeric(3)
for (gi in seq_along(cc)){
  obj <- cc[[gi]]; idents <- levels(obj@idents)
  if (!"NGL" %in% obj@netP$pathways) next
  pi <- which(obj@netP$pathways == "NGL")
  ai <- which(idents %in% ASTROS); ni <- which(idents == "Excitatory_Neuron")
  ngl_sum[gi] <- sum(obj@netP$prob[ai, ni, pi], na.rm=TRUE)
}
ngl_pct <- ngl_sum / ngl_sum[1] * 100
cat("NGL 总强度:", round(ngl_sum, 4), "  % vs Sham:", round(ngl_pct,1), "\n")

# ===== 数据 2 (Panel B): NRXN 多度量 =====
# B-上: NRXN LR-pair 平均强度 (Astro→Neuron)
nrxn_mean <- numeric(3)
for (gi in seq_along(cc)){
  obj <- cc[[gi]]
  if (!"NRXN" %in% obj@netP$pathways) next
  net <- subsetCommunication(obj, signaling="NRXN")
  sub <- net[net$source %in% ASTROS & net$target %in% NEURONS, ]
  nrxn_mean[gi] <- mean(sub$prob, na.rm=TRUE)
}
nrxn_mean_pct <- nrxn_mean / nrxn_mean[1] * 100
cat("NRXN mean prob:", round(nrxn_mean, 4), "  % vs Sham:", round(nrxn_mean_pct,1), "\n")

# B-下: Astro_C1 NRXN 总输出强度 (sum outgoing)
nrxn_c1 <- numeric(3)
for (gi in seq_along(cc)){
  obj <- cc[[gi]]; idents <- levels(obj@idents)
  if (!"NRXN" %in% obj@netP$pathways) next
  pi <- which(obj@netP$pathways == "NRXN")
  ai <- which(idents == "Astro_C1")
  nrxn_c1[gi] <- sum(obj@netP$prob[ai, , pi], na.rm=TRUE)
}
nrxn_c1_pct <- nrxn_c1 / nrxn_c1[1] * 100
cat("Astro_C1 NRXN outgoing:", round(nrxn_c1, 4), "  % vs Sham:", round(nrxn_c1_pct,1), "\n")

# ===== 数据 3 (Panel C): 亚群 × 通路 log2FC heatmap =====
PATHS <- c("NRXN","NCAM","NGL")
B_data <- array(0, dim=c(length(ASTROS), length(PATHS), 3),
                dimnames=list(ASTROS, PATHS, c("Sham","CCI","CCI+KO")))
get_prob <- function(obj, a, n, pp){
  idents <- levels(obj@idents)
  if (! pp %in% obj@netP$pathways) return(0)
  if (! a %in% idents || ! n %in% idents) return(0)
  pi <- which(obj@netP$pathways == pp)
  ai <- which(idents==a); ni <- which(idents==n)
  obj@netP$prob[ai, ni, pi]
}
for (ai in seq_along(ASTROS)) for (pi in seq_along(PATHS)) for (gi in seq_along(cc))
  B_data[ai, pi, gi] <- get_prob(cc[[gi]], ASTROS[ai], "Excitatory_Neuron", PATHS[pi])
pseudo <- 0.001
log2fc <- log2((B_data[,,c("CCI","CCI+KO")] + pseudo) /
               array(B_data[,,"Sham"] + pseudo, dim=c(3,3,2)))

# ===== 绘图 =====
png(OUT, width=3600, height=1400, res=300, type="cairo", family="Microsoft YaHei")
layout(matrix(c(1,2,3), 1, 3), widths=c(1, 1.1, 1.3))
par(mar=c(5.5, 5.5, 4.2, 1.0), mgp=c(3.0, 0.7, 0))

cols3 <- c(navy, accent, teal)
lbl3 <- c("Sham","CCI","CCI+TGF-βRI-KO")

# --- Panel A: NGL ---
bp <- barplot(ngl_pct, col=cols3, border=NA, names.arg=lbl3,
              ylim=c(0, max(ngl_pct)*1.20),
              ylab="NGL 通讯强度 (% vs Sham)", cex.names=0.85, font.axis=2)
abline(h=100, col=gray1, lty=2)
for (i in 1:3)
  text(bp[i], ngl_pct[i]+max(ngl_pct)*0.04, sprintf("%.0f%%", ngl_pct[i]),
       cex=1.0, font=2, col=cols3[i])
mtext("A", side=3, line=2.3, adj=-0.30, font=2, cex=2.0, col="black")
mtext("NGL 通路 Astro→Exc Neuron (强支持)", side=3, line=0.4, cex=0.92, col=ink)

# --- Panel B: NRXN 多度量 (2 个柱图叠加) ---
par(mar=c(5.5, 5.5, 4.2, 1.0))
# 矩阵: 3 行 (groups) × 2 列 (metrics) → beside=TRUE 每 metric 一组 3 bars
b_data <- cbind("LR-pair 平均强度"=c(Sham=nrxn_mean_pct[1], CCI=nrxn_mean_pct[2], KO=nrxn_mean_pct[3]),
                "Astro_C1 输出强度" =c(Sham=nrxn_c1_pct[1],   CCI=nrxn_c1_pct[2],   KO=nrxn_c1_pct[3]))
bp2 <- barplot(b_data, beside=TRUE, col=cols3, border=NA,
               names.arg=rep("", ncol(b_data)),
               ylim=c(0, max(b_data)*1.25), space=c(0.2, 0.8),
               ylab="NRXN 通讯 (% vs Sham)", font.axis=2)
abline(h=100, col=gray1, lty=2)
metric_centers <- colMeans(bp2)
axis(1, at=metric_centers, labels=colnames(b_data), tick=FALSE,
     font=2, cex.axis=0.92, line=-0.3)
for (j in 1:2) for (i in 1:3)
  text(bp2[i,j], b_data[i,j]+max(b_data)*0.04, sprintf("%.0f%%", b_data[i,j]),
       cex=0.85, font=2, col=cols3[i])
legend("topleft", legend=lbl3, fill=cols3, bty="o", box.col="white",
       bg=adjustcolor("white",0.9), cex=0.78, inset=c(0.02, 0.02))
mtext("B", side=3, line=2.3, adj=-0.20, font=2, cex=2.0, col="black")
mtext("NRXN 通路多度量验证 (中等支持)", side=3, line=0.4, cex=0.92, col=ink)

# --- Panel C: 亚群×通路 log2FC heatmap ---
par(mar=c(5.5, 8.0, 4.2, 3.5))
hm <- matrix(NA, length(ASTROS)*length(PATHS), 2)
row_lbls <- character(length(ASTROS)*length(PATHS))
r <- 1
for (a in ASTROS) for (p in PATHS){
  row_lbls[r] <- paste0(a, "→", p)
  hm[r, 1] <- log2fc[a, p, "CCI"]
  hm[r, 2] <- log2fc[a, p, "CCI+KO"]
  r <- r + 1
}
rownames(hm) <- row_lbls
colnames(hm) <- c("CCI/Sham","CCI+KO/Sham")
nr <- nrow(hm); nc <- ncol(hm)
zlim <- max(abs(hm), na.rm=TRUE)
brk <- seq(-zlim, zlim, length.out=101)
cols_hm <- colorRampPalette(c("#2166AC","#67A9CF","#F7F7F7","#EF8A62","#B2182B"))(100)
image(1:nc, 1:nr, t(hm)[, nr:1], col=cols_hm, breaks=brk, axes=FALSE, xlab="", ylab="")
rect(0.5, 0.5, nc+0.5, nr+0.5, border="black", lwd=1.5)
axis(1, at=1:nc, labels=colnames(hm), tick=FALSE, font=2, cex.axis=0.95, line=-0.3)
for (i in 1:nr)
  axis(2, at=i, labels=rev(row_lbls)[i], las=2, tick=FALSE, cex.axis=0.85, font=2)
hm_t <- t(hm)[, nr:1]
for (j in 1:nc) for (i in 1:nr){
  v <- hm_t[j, i]; if (is.na(v)) next
  txt_col <- if (abs(v) > zlim*0.55) "white" else ink
  text(j, i, sprintf("%+.1f", v), cex=0.85, col=txt_col, font=2)
}
for (k in 1:2) abline(h=k*length(PATHS)+0.5, col="black", lwd=1.2, lty=2)
par(xpd=NA)
cb_x <- nc + 0.85; cb_w <- 0.20
cb_y <- seq(0.5, nr+0.5, length.out=length(cols_hm)+1)
for (k in 1:length(cols_hm))
  rect(cb_x, cb_y[k], cb_x+cb_w, cb_y[k+1], col=cols_hm[k], border=NA)
text(cb_x+cb_w/2, nr+0.85, "log2FC", cex=0.78, font=2)
text(cb_x+cb_w+0.05, c(0.5, (nr+1)/2, nr+0.5),
     c(sprintf("%+.1f",-zlim), "0", sprintf("%+.1f",zlim)),
     cex=0.78, adj=0)
par(xpd=FALSE)
mtext("C", side=3, line=2.3, adj=-0.20, font=2, cex=2.0, col="black")
mtext("亚群层次异质性 (log2FC, Astro_sub→Exc Neuron)", side=3, line=0.4, cex=0.92, col=ink)

dev.off()
cat("\n图_selftest_v3 已生成 ->", OUT, "\n")
