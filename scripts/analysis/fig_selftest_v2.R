#!/usr/bin/env Rscript
# 研究基础新增图: 自测 snRNA 3 组 (Sham/CCI/CCI+TGFβRI-KO)
# Astro 亚群 → Neuron 黏附通讯
# A: 通路总强度 3 组对比 (NRXN/NCAM/NGL)
# B: 亚群 × 通路热图 (% vs Sham), 直接展示亚群异质性
suppressMessages({library(CellChat)})
source("D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/style_lib.R")

D <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/wet_lab_data/results/cellchat_nrxn"
OUT <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/fig_selftest_v2.png"

cc <- list(
  Sham = readRDS(file.path(D,"cellchat_nrxn_Sham.rds")),
  CCI  = readRDS(file.path(D,"cellchat_nrxn_CCI_Scr.rds")),
  KO   = readRDS(file.path(D,"cellchat_nrxn_CCI_TGFbetaR1_KO.rds"))
)

ASTROS <- c("Astro_C0","Astro_C1","Astro_C2")
NEURONS <- c("Excitatory_Neuron","Inhibitory_Neuron")
PATHS  <- c("NRXN","NCAM","NGL")  # NEGR 在所有亚群均=0, 不展示

# 提取 Astro_sub × Neuron × Pathway × Group 强度
get_prob <- function(obj, a, n, pp){
  idents <- levels(obj@idents)
  if (! pp %in% obj@netP$pathways) return(0)
  if (! a %in% idents || ! n %in% idents) return(0)
  pi <- which(obj@netP$pathways == pp)
  ai <- which(idents==a); ni <- which(idents==n)
  obj@netP$prob[ai, ni, pi]
}

# === Panel A 数据: 各通路 Astro→Exc_Neuron 总强度 (sum over astro_sub) 3 组 ===
A_data <- matrix(0, length(PATHS), 3, dimnames=list(PATHS, c("Sham","CCI","CCI+KO")))
for (pi in seq_along(PATHS)){
  for (gi in seq_along(cc)){
    s <- 0
    for (a in ASTROS) s <- s + get_prob(cc[[gi]], a, "Excitatory_Neuron", PATHS[pi])
    A_data[pi, gi] <- s
  }
}
cat("=== A 数据 (Astro→Exc_Neuron 总强度) ===\n"); print(round(A_data, 4))
A_pct <- sweep(A_data, 1, A_data[,"Sham"], "/") * 100
A_pct[is.na(A_pct) | is.infinite(A_pct)] <- 0
cat("\n=== A 百分比 (% vs Sham) ===\n"); print(round(A_pct, 1))

# === Panel B 数据: Astro_sub × pathway × group (% vs Sham 内同亚群同通路) ===
B_data <- array(0, dim=c(length(ASTROS), length(PATHS), 3),
                dimnames=list(ASTROS, PATHS, c("Sham","CCI","CCI+KO")))
for (ai in seq_along(ASTROS)){
  for (pi in seq_along(PATHS)){
    for (gi in seq_along(cc)){
      B_data[ai, pi, gi] <- get_prob(cc[[gi]], ASTROS[ai], "Excitatory_Neuron", PATHS[pi])
    }
  }
}
cat("\n=== B 原始数据 ===\n"); print(round(B_data, 4))

# 计算 log2(FC vs Sham); 0 → -Inf 处理: 加一个小 pseudocount
pseudo <- 0.001
B_l2fc <- log2((B_data[,,c("CCI","CCI+KO")] + pseudo) /
               array(B_data[,,"Sham"] + pseudo, dim=c(3,3,2)))
cat("\n=== B log2FC vs Sham ===\n"); print(round(B_l2fc, 2))

# ===== 绘图 =====
png(OUT, width=3200, height=1500, res=300, type="cairo", family="Microsoft YaHei")
par(mfrow=c(1,2), mar=c(5.5, 5.5, 3.5, 1.0), mgp=c(3.0, 0.7, 0))

# --- Panel A: 分组柱状 ---
cols3 <- c(navy, accent, teal)
bp <- barplot(t(A_pct), beside=TRUE, col=cols3, border=NA,
              ylim=c(0, max(A_pct, na.rm=TRUE)*1.25), space=c(0.2, 0.8),
              ylab="Astro→Exc Neuron 通讯强度 (% vs Sham)",
              names.arg=rep("", length(PATHS)))
abline(h=100, col=gray1, lty=2, lwd=1.2)
axis(1, at=colMeans(bp), labels=PATHS, tick=FALSE, font=2, cex.axis=1.05, line=-0.2)
axis(2, cex.axis=0.95)
for (j in 1:length(PATHS)){
  for (i in 1:3){
    text(bp[i,j], A_pct[j,i]+max(A_pct, na.rm=TRUE)*0.04,
         sprintf("%.0f%%", A_pct[j,i]), cex=0.78, font=2, col=cols3[i])
  }
}
legend("topright", legend=c("Sham","CCI","CCI+TGF-βRI-KO"),
       fill=cols3, bty="o", box.col="white", bg=adjustcolor("white",0.92),
       cex=0.85, inset=c(0.01, 0.01))
mtext("A", side=3, line=1.4, adj=-0.18, font=2, cex=2.0, col="black")
mtext("Astro→Exc Neuron 主要黏附通路时序变化", side=3, line=0.3, cex=0.92, col=ink)

# --- Panel B: 亚群异质性 log2FC heatmap (CCI vs Sham + KO vs Sham) ---
par(mar=c(5.5, 8.0, 3.5, 4.0))
# 拼成 (Astro × Path) × (CCI, KO) 矩阵
hm <- matrix(NA, length(ASTROS)*length(PATHS), 2)
row_lbls <- character(length(ASTROS)*length(PATHS))
r <- 1
for (a in ASTROS){
  for (p in PATHS){
    row_lbls[r] <- paste0(a, "→", p)
    hm[r, 1] <- B_l2fc[a, p, "CCI"]
    hm[r, 2] <- B_l2fc[a, p, "CCI+KO"]
    r <- r + 1
  }
}
rownames(hm) <- row_lbls
colnames(hm) <- c("CCI/Sham","CCI+KO/Sham")

nr <- nrow(hm); nc <- ncol(hm)
zlim <- max(abs(hm), na.rm=TRUE)
brk <- seq(-zlim, zlim, length.out=101)
cols_hm <- colorRampPalette(c("#2166AC","#67A9CF","#F7F7F7","#EF8A62","#B2182B"))(100)
image(1:nc, 1:nr, t(hm)[, nr:1], col=cols_hm, breaks=brk, axes=FALSE,
      xlab="", ylab="", useRaster=FALSE)
rect(0.5, 0.5, nc+0.5, nr+0.5, border="black", lwd=1.5)
# X 轴
axis(1, at=1:nc, labels=colnames(hm), tick=FALSE, font=2, cex.axis=1.0, line=-0.2)
# Y 轴
for (i in 1:nr){
  lab <- rev(row_lbls)[i]
  axis(2, at=i, labels=lab, las=2, tick=FALSE, cex.axis=0.92, font=2)
}
# 数值标注
hm_t <- t(hm)[, nr:1]
for (j in 1:nc){
  for (i in 1:nr){
    v <- hm_t[j, i]
    if (is.na(v)) next
    txt_col <- if (abs(v) > zlim*0.55) "white" else ink
    text(j, i, sprintf("%+.1f", v), cex=0.92, col=txt_col, font=2)
  }
}
# 亚群分隔线 (每 3 行一组)
for (k in 1:2) abline(h=k*length(PATHS)+0.5, col="black", lwd=1.2, lty=2)
# color bar
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
mtext("B", side=3, line=1.4, adj=-0.20, font=2, cex=2.0, col="black")
mtext("亚群层次异质性 (log2 倍数变化)", side=3, line=0.3, cex=0.92, col=ink)

dev.off()
cat("\n图_selftest_v2 已生成 ->", OUT, "\n")
