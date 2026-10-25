#!/usr/bin/env Rscript
# 研究基础新增图: 自测 snRNA-seq 三组 (Sham/CCI_Scr/CCI_TGFbetaR1_KO)
# Astro → Neuron 的 NRXN/NCAM/NEGR/NGL 4 大黏附通讯通路强度对比
# 支持立论依据中"45-55%下降 + KO 回升"的声明
suppressMessages({library(CellChat); library(Matrix)})
source("D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/style_lib.R")

D <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/wet_lab_data/results/cellchat_nrxn"
OUT <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/fig_selftest_adhesion.png"

cc <- list(
  Sham = readRDS(file.path(D,"cellchat_nrxn_Sham.rds")),
  CCI  = readRDS(file.path(D,"cellchat_nrxn_CCI_Scr.rds")),
  KO   = readRDS(file.path(D,"cellchat_nrxn_CCI_TGFbetaR1_KO.rds"))
)
cat("loaded 3 CellChat objects\n")

PATHWAYS <- c("NRXN","NCAM","NEGR","NGL")

# 提取每组每通路 Astro → Neuron 的通讯强度
# 步骤: 用 cc@netP$prob (array: source × target × pathway)
# Astro = 所有含"Astro"或"OPC"的 source idents
# Neuron = 所有含"CA_/CA1/CA3/DG_granule/Inh_/Camk2a/Slc17a6/PV/SST/VIP/NPY/Lamp5/Reactive"等 target

extract_path_strength <- function(obj, pathway, src_pattern, tgt_pattern){
  idents <- levels(obj@idents)
  src_idx <- grep(src_pattern, idents, ignore.case=TRUE)
  tgt_idx <- grep(tgt_pattern, idents, ignore.case=TRUE)
  if (! pathway %in% obj@netP$pathways) return(NA)
  pp <- which(obj@netP$pathways == pathway)
  prob_mat <- obj@netP$prob[,,pp]  # src × tgt
  v <- prob_mat[src_idx, tgt_idx]
  sum(v, na.rm=TRUE)   # 总通讯强度
}

# 模式: Astro = "Astro|Reactive|DAM"; Neuron = 各 Neuron 类型
# 检查 idents
cat("\n=== idents per group ===\n")
for (n in names(cc)) {
  cat(n, ":\n")
  print(head(levels(cc[[n]]@idents), 30))
}

# cellchat_nrxn 系列 idents: Astro_C0/C1/C2 + Excitatory_Neuron + Inhibitory_Neuron
SRC <- "^Astro_C"            # 3 个 Astro 亚群
TGT <- "Neuron"               # Excitatory_Neuron + Inhibitory_Neuron

# 提取
M <- matrix(NA, length(PATHWAYS), 3, dimnames=list(PATHWAYS, c("Sham","CCI","KO")))
for (gi in seq_along(cc)){
  g <- names(cc)[gi]
  for (pi in seq_along(PATHWAYS)){
    p <- PATHWAYS[pi]
    M[pi, gi] <- extract_path_strength(cc[[g]], p, SRC, TGT)
  }
}
cat("\n=== Pathway strength (sum prob, Astro→Neuron) ===\n")
print(round(M, 4))

# 归一化到 Sham=100%
pct <- t(t(M) / M[,"Sham"] * 100)   # 错误, 维度需要 colMeans... 直接除
pct <- sweep(M, 1, M[,"Sham"], "/") * 100
cat("\n=== %  (vs Sham=100%) ===\n")
print(round(pct, 1))

# Save data csv
write.csv(M,   "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/selftest_adhesion_strength.csv")
write.csv(pct, "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/selftest_adhesion_pct.csv")

# ===== 绘图: 2 面板 =====
# A: 4 通路 × 3 组 分组柱状 (pct vs Sham)
# B: NRXN 通路 主导 LR pairs (NRXN1-NLGN1 等) 在 3 组的概率热图

png(OUT, width=3200, height=1400, res=300, type="cairo", family="Microsoft YaHei")
par(mfrow=c(1,2), mar=c(6.0, 5.5, 3.5, 1.0), mgp=c(3.0, 0.7, 0))

# --- Panel A: 分组柱状 ---
n_path <- nrow(pct); n_grp <- ncol(pct)
cols3 <- c(navy, accent, teal)
bp <- barplot(t(pct), beside=TRUE, col=cols3, border=NA,
              ylim=c(0, max(pct)*1.20), space=c(0.2, 0.8),
              ylab="Astro→Neuron 通讯强度 (% vs Sham)",
              names.arg=rep("", n_path))
abline(h=100, col=gray1, lty=2, lwd=1.2)
# 各通路名 (位于 group 组下方)
axis(1, at=colMeans(bp), labels=PATHWAYS, tick=FALSE, font=2, cex.axis=1.0, line=-0.2)
axis(2, cex.axis=0.95)
# 数值标注 on top of each bar
for (j in 1:n_path){
  for (i in 1:n_grp){
    text(bp[i,j], pct[j,i]+max(pct)*0.04, sprintf("%.0f%%", pct[j,i]),
         cex=0.78, font=2, col=cols3[i])
  }
}
# 图例
legend("topright", legend=c("Sham","CCI","CCI+TGF-βRI-KO"),
       fill=cols3, bty="o", box.col="white", bg=adjustcolor("white",0.9),
       cex=0.88, inset=c(0.01, 0.01))
mtext("A", side=3, line=1.4, adj=-0.18, font=2, cex=2.0, col="black")
mtext("Astro→Neuron 突触黏附通路强度 (4 大通路 × 3 组)",
      side=3, line=0.4, cex=0.92, col=ink)

# --- Panel B: NRXN 主导 LR 对 概率 (列表) ---
# 从 prob array 中取 NRXN 通路, 列出 top LR pairs
get_top_LR <- function(obj, pathway, src_pattern, tgt_pattern, top=6){
  idents <- levels(obj@idents)
  src_idx <- grep(src_pattern, idents, ignore.case=TRUE)
  tgt_idx <- grep(tgt_pattern, idents, ignore.case=TRUE)
  # LR pair 在 net (不是 netP) 里
  net <- subsetCommunication(obj, signaling=pathway)
  net <- net[net$source %in% idents[src_idx] & net$target %in% idents[tgt_idx], ]
  agg <- aggregate(prob ~ interaction_name_2, net, sum)
  agg <- agg[order(-agg$prob), ]
  head(agg, top)
}
# 取 Sham 的 top NRXN LR pairs 作为列名
sham_lr <- get_top_LR(cc$Sham, "NRXN", SRC, TGT, top=6)
lr_names <- sham_lr$interaction_name_2
# 提取 3 组每组对这些 LR 的总 prob
lr_M <- matrix(NA, length(lr_names), 3, dimnames=list(lr_names, c("Sham","CCI","KO")))
for (gi in seq_along(cc)){
  g <- names(cc)[gi]
  net <- subsetCommunication(cc[[g]], signaling="NRXN")
  idents <- levels(cc[[g]]@idents)
  src_idx <- grep(SRC, idents, ignore.case=TRUE)
  tgt_idx <- grep(TGT, idents, ignore.case=TRUE)
  net <- net[net$source %in% idents[src_idx] & net$target %in% idents[tgt_idx], ]
  for (li in seq_along(lr_names)){
    ln <- lr_names[li]
    lr_M[li, gi] <- sum(net$prob[net$interaction_name_2 == ln], na.rm=TRUE)
  }
}
cat("\n=== NRXN top LR pairs (sum prob) ===\n"); print(round(lr_M, 4))

# heatmap
zlim <- max(abs(lr_M), na.rm=TRUE)
brk <- seq(0, zlim, length.out=101)
cols_hm <- colorRampPalette(c("#F7F7F7","#F4A582","#D6604D","#B2182B"))(100)
par(mar=c(6.0, 9.0, 3.5, 3.5))
nr <- nrow(lr_M); nc <- ncol(lr_M)
image(1:nc, 1:nr, t(lr_M)[, nr:1], col=cols_hm, breaks=brk, axes=FALSE,
      xlab="", ylab="", useRaster=FALSE)
rect(0.5, 0.5, nc+0.5, nr+0.5, border="black", lwd=1.5)
axis(1, at=1:nc, labels=c("Sham","CCI","CCI+KO"), tick=FALSE,
     font=2, cex.axis=1.0, line=-0.2)
for (i in 1:nr){
  ln <- rev(lr_names)[i]
  axis(2, at=i, labels=ln, las=2, tick=FALSE, cex.axis=0.85, font=2)
}
# 数值标注
for (j in 1:nc){
  for (i in 1:nr){
    val <- t(lr_M)[, nr:1][j, i]
    txt_col <- if (val > zlim*0.5) "white" else ink
    text(j, i, sprintf("%.3f", val), cex=0.80, col=txt_col, font=2)
  }
}
# color bar
par(xpd=NA)
cb_x <- nc + 0.85; cb_w <- 0.18
cb_y <- seq(0.5, nr+0.5, length.out=length(cols_hm)+1)
for (k in 1:length(cols_hm))
  rect(cb_x, cb_y[k], cb_x+cb_w, cb_y[k+1], col=cols_hm[k], border=NA)
text(cb_x+cb_w/2, nr+0.85, "Σ prob", cex=0.78, font=2)
text(cb_x+cb_w+0.05, c(0.5, nr+0.5),
     c("0", sprintf("%.2f", zlim)), cex=0.78, adj=0)
par(xpd=FALSE)

mtext("B", side=3, line=1.4, adj=-0.20, font=2, cex=2.0, col="black")
mtext("NRXN 通路主导配体-受体对 强度热图 (Astro→Neuron)",
      side=3, line=0.4, cex=0.92, col=ink)

dev.off()
cat("\n图_selftest 已生成 ->", OUT, "\n")
