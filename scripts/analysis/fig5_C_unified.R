#!/usr/bin/env Rscript
# 图5 C 重做 fEPSP Slope (%) - 4 组 bar (无散点)
source("D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/style_lib.R")

OUT <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/fig5_C_unified.png"

labels4 <- c("Sham","TBI","TBI+Ab(20pg/mg)","TBI+Ab(60pg/mg)")
means <- c(155, 132, 142, 145)
errs  <- c(2.5, 3, 4, 3.5)

sig_brk <- list(
  list(from=1, to=2, p=0.005),
  list(from=2, to=3, p=0.03),
  list(from=2, to=4, p=0.03)
)

png(OUT, width=1500, height=1500, res=300, type="cairo", family="Microsoft YaHei")
par(mar=c(7.5, 5.0, 3.2, 1.2), mgp=c(2.8, 0.7, 0))

cols4 <- c(navy, accent, steel, teal)
unified_bar(means, errs, labels4, cols4, dots=NULL,
            sig_brackets=sig_brk, ylab="fEPSP Slope (%)",
            ylim=c(0, 200), h0=FALSE)
mtext("C", side=3, line=1.4, adj=-0.20, font=2, cex=2.2, col="black")

dev.off()
cat("图5 C unified 已生成\n")
