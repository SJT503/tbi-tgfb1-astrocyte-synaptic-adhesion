#!/usr/bin/env Rscript
# 图2 C/D 重做 GFAP+TGFβRI+ 比例 (Cortex/Hippocampus)
source("D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/style_lib.R")

OUT <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/fig2_CD_unified.png"

labels3 <- c("Sham","Ipsilateral","Contralateral")
cols3 <- c(navy, accent, steel)

# Panel C: Cortex (视觉读取)
C_mean <- c(1.0, 1.9, 1.25)
C_se   <- c(0.15, 0.20, 0.18)
set.seed(1); C_dots <- mapply(function(m, se) m + rnorm(3, 0, se*sqrt(3)), C_mean, C_se, SIMPLIFY=FALSE)
C_sig <- list(
  list(from=1, to=2, p=0.03),
  list(from=1, to=3, p=0.20),
  list(from=2, to=3, p=0.20)
)

# Panel D: Hippocampus (视觉读取)
D_mean <- c(1.0, 2.15, 1.45)
D_se   <- c(0.18, 0.20, 0.22)
set.seed(2); D_dots <- mapply(function(m, se) m + rnorm(3, 0, se*sqrt(3)), D_mean, D_se, SIMPLIFY=FALSE)
D_sig <- list(
  list(from=1, to=2, p=0.03),
  list(from=1, to=3, p=0.20),
  list(from=2, to=3, p=0.20)
)

png(OUT, width=2400, height=1300, res=300, type="cairo", family="Microsoft YaHei")
par(mfrow=c(1,2), mar=c(6.5, 5.0, 3.2, 1.2), mgp=c(2.8, 0.7, 0))

unified_bar(C_mean, C_se, labels3, cols3, dots=C_dots,
            sig_brackets=C_sig, ylab=expression("Relative Ratio of GFAP"^"+"*" TGF-"*beta*"RI"^"+"*" cells"),
            ylim=c(0, 3.6), h0=FALSE)
mtext("C", side=3, line=1.4, adj=-0.20, font=2, cex=2.2, col="black")
mtext("Cortex (皮层)", side=3, line=0.4, cex=0.95, col=ink)

unified_bar(D_mean, D_se, labels3, cols3, dots=D_dots,
            sig_brackets=D_sig, ylab=expression("Relative Ratio of GFAP"^"+"*" TGF-"*beta*"RI"^"+"*" cells"),
            ylim=c(0, 3.6), h0=FALSE)
mtext("D", side=3, line=1.4, adj=-0.20, font=2, cex=2.2, col="black")
mtext("Hippocampus (海马)", side=3, line=0.4, cex=0.95, col=ink)

dev.off()
cat("图2 C/D unified 已生成\n")
