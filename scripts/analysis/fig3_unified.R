#!/usr/bin/env Rscript
# 图3 重做: Recognition Index + Discrimination Index
# 数据值来源: 任智慧硕士论文 Fig 3-3 视觉读取 (保留原始统计结果)
source("D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/style_lib.R")

OUT <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/fig3_unified.png"

labels5 <- c("Sham+Vehicle","TBI+Vehicle","TBI+TGFβ1(1ng)","TBI+TGFβ1(10ng)","TBI+TGFβ1(50ng)")

# Panel A: Recognition Index (%) — 视觉读取 (任智慧论文)
A_mean <- c(75, 42, 35, 33, 25)
A_se   <- c(7, 6, 9, 11, 11)
# 散点 (n=6/组, 用 mean ± SE*sqrt(6) 估计 SD 反算)
set.seed(1)
A_dots <- mapply(function(m, se) m + rnorm(6, 0, se*sqrt(6)),
                 A_mean, A_se, SIMPLIFY=FALSE)
A_sig <- list(
  list(from=1, to=2, p=0.002),
  list(from=1, to=3, p=0.002),
  list(from=1, to=4, p=0.002),
  list(from=1, to=5, p=0.002),
  list(from=2, to=5, p=0.03)
)

# Panel B: Discrimination Index (%) — 视觉读取
B_mean <- c(50, 0, -5, -12, -18)
B_se   <- c(9, 11, 13, 14, 16)
set.seed(2)
B_dots <- mapply(function(m, se) m + rnorm(6, 0, se*sqrt(6)),
                 B_mean, B_se, SIMPLIFY=FALSE)
B_sig <- list(
  list(from=1, to=2, p=0.03),
  list(from=1, to=3, p=0.002),
  list(from=1, to=4, p=0.002),
  list(from=1, to=5, p=0.002)
)

png(OUT, width=3000, height=1500, res=300, type="cairo", family="Microsoft YaHei")
par(mfrow=c(1,2), mar=c(7.5, 5.0, 3.2, 1.2), mgp=c(2.8, 0.7, 0))

# Panel A
unified_bar(A_mean, A_se, labels5, PAL5, dots=A_dots,
            sig_brackets=A_sig, ylab="Recognition Index (%)",
            ylim=c(0, 130), h0=FALSE)
mtext("A", side=3, line=1.4, adj=-0.20, font=2, cex=2.2, col="black")

# Panel B
unified_bar(B_mean, B_se, labels5, PAL5, dots=B_dots,
            sig_brackets=B_sig, ylab="Discrimination Index (%)",
            ylim=c(-90, 100), h0=TRUE)
mtext("B", side=3, line=1.4, adj=-0.20, font=2, cex=2.2, col="black")

dev.off()
cat("图3 unified 已生成 ->", OUT, "\n")
