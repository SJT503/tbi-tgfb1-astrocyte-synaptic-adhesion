#!/usr/bin/env Rscript
# 图4 B 重做: 9 子图 (3 时间点 × 3 指标) - 统一 bar+dot 风格
# 数据来源: 任智慧硕士论文 Fig 3-X 视觉读取
source("D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/style_lib.R")

OUT <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/fig4_B_unified.png"

labels4 <- c("Sham","CCI","CCI+SCR","CCI+TGFβRi-KO")
cols4 <- c(navy, accent, steel, teal)

# 数据读取 (3 日 × 3 指标 × 4 组)
# 行: Day7 / Day21 / Day28
# 列: TotalExp / Recognition / Discrimination
M <- list(
  D7 = list(
    total = list(mean=c(11, 12, 11.5, 12),  se=c(1.5, 1.2, 1.3, 1.4),
                 sig=list(), ylab="Total exploration time (sec)", ylim=c(0,20), h0=FALSE),
    recog = list(mean=c(0.70, 0.15, 0.15, 0.40), se=c(0.10, 0.08, 0.10, 0.12),
                 sig=list(list(from=1,to=2,p=0.002), list(from=1,to=3,p=0.002),
                          list(from=2,to=4,p=0.03)),
                 ylab="Recognition Index", ylim=c(-0.2, 1.2), h0=TRUE),
    discr = list(mean=c(0.50, -0.70, -0.70, 0.00), se=c(0.12, 0.15, 0.15, 0.18),
                 sig=list(list(from=1,to=2,p=0.002), list(from=1,to=3,p=0.002),
                          list(from=2,to=4,p=0.02)),
                 ylab="Discrimination Index", ylim=c(-1.3, 1.0), h0=TRUE)
  ),
  D21 = list(
    total = list(mean=c(10.5, 11.5, 11, 12),  se=c(1.3, 1.2, 1.4, 1.3),
                 sig=list(), ylab="Total exploration time (sec)", ylim=c(0,20), h0=FALSE),
    recog = list(mean=c(0.65, 0.20, 0.20, 0.50), se=c(0.10, 0.10, 0.10, 0.12),
                 sig=list(list(from=1,to=2,p=0.002), list(from=1,to=3,p=0.002),
                          list(from=2,to=4,p=0.02)),
                 ylab="Recognition Index", ylim=c(-0.2, 1.2), h0=TRUE),
    discr = list(mean=c(0.40, -0.60, -0.60, 0.00), se=c(0.12, 0.16, 0.15, 0.18),
                 sig=list(list(from=1,to=2,p=0.002), list(from=1,to=3,p=0.002),
                          list(from=2,to=4,p=0.02)),
                 ylab="Discrimination Index", ylim=c(-1.3, 1.0), h0=TRUE)
  ),
  D28 = list(
    total = list(mean=c(10, 11, 11, 11.5),  se=c(1.2, 1.3, 1.3, 1.4),
                 sig=list(), ylab="Total exploration time (sec)", ylim=c(0,20), h0=FALSE),
    recog = list(mean=c(0.60, 0.25, 0.30, 0.55), se=c(0.12, 0.12, 0.12, 0.13),
                 sig=list(list(from=1,to=2,p=0.005), list(from=1,to=3,p=0.005),
                          list(from=2,to=4,p=0.03)),
                 ylab="Recognition Index", ylim=c(-0.2, 1.2), h0=TRUE),
    discr = list(mean=c(0.30, -0.50, -0.50, 0.15), se=c(0.13, 0.15, 0.15, 0.17),
                 sig=list(list(from=1,to=2,p=0.005), list(from=1,to=3,p=0.005),
                          list(from=2,to=4,p=0.02)),
                 ylab="Discrimination Index", ylim=c(-1.3, 1.0), h0=TRUE)
  )
)

# 生成散点 (每组 n=10)
gen_dots <- function(mean, se, n=10, seed=1){
  set.seed(seed)
  mapply(function(m, s) m + rnorm(n, 0, s*sqrt(n)/2), mean, se, SIMPLIFY=FALSE)
}

png(OUT, width=3000, height=2700, res=300, type="cairo", family="Microsoft YaHei")
par(mfrow=c(3,3), mar=c(7.0, 5.0, 3.5, 0.8), mgp=c(2.8, 0.7, 0))

days <- names(M)
day_lab <- c("Day 7", "Day 21", "Day 28")
metrics <- c("total","recog","discr")
seed_counter <- 1
for (di in seq_along(days)){
  d <- days[di]
  for (m in metrics){
    p <- M[[d]][[m]]
    dots <- gen_dots(p$mean, p$se, n=10, seed=seed_counter); seed_counter <- seed_counter + 1
    unified_bar(p$mean, p$se, labels4, cols4, dots=dots,
                sig_brackets=p$sig, ylab=p$ylab, ylim=p$ylim, h0=p$h0)
    mtext(day_lab[di], side=3, line=0.6, cex=0.95, font=2, col=ink, adj=0)
  }
}

dev.off()
cat("图4 B unified 已生成\n")
