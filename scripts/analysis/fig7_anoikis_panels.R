#!/usr/bin/env Rscript
# 图7 panel C/D: 神经元 anoikis 通路时序
A <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/anoikis_out"
OUT <- A
agg <- read.csv(file.path(A,"anoikis_per_sample.csv"))
agg$timepoint <- factor(agg$timepoint, levels=c("NoTBI","D1","D7","D28"))

navy<-"#1a3a5c"; mid<-"#3f6896"; steel<-"#7795b5"; accent<-"#b5563f"; teal<-"#2a6a5e"
ink<-"#2a3a4a"

draw <- function() {
  par(mfrow=c(1,2), mar=c(4.6,5.0,3.4,1.0), mgp=c(2.7,0.7,0),
      family="Microsoft YaHei", cex.axis=0.94, cex.lab=1.05)

  ## Panel C: 凋亡执行通路 死亡评分 时序
  m  <- tapply(agg$death, agg$timepoint, mean)
  se <- tapply(agg$death, agg$timepoint, function(x) sd(x)/sqrt(length(x)))
  bp <- barplot(m, col=c(navy,accent,accent,accent), border=NA,
                ylim=c(min(m-se)*1.6, max(m+se)*2.2),
                ylab="凋亡执行通路评分", main="")
  abline(h=0, col="#9aa7b3", lty=2)
  arrows(bp, m-se, bp, m+se, angle=90, code=3, length=0.05, col=ink, lwd=1.4)
  set.seed(1)
  for (i in seq_along(levels(agg$timepoint))) {
    tp <- levels(agg$timepoint)[i]; y <- agg$death[agg$timepoint==tp]
    points(jitter(rep(bp[i],length(y)),amount=0.10), y, pch=21, bg="white", col=ink, cex=0.95)
  }
  text(bp[2], m[2]+se[2]+0.025, "*", cex=1.9, col=accent, font=2)
  text(bp[3], m[3]+se[3]+0.025, "*", cex=1.9, col=accent, font=2)
  text(bp[4], m[4]+se[4]+0.025, "*", cex=1.9, col=accent, font=2)
  mtext("C", side=3, line=1.6, adj=-0.20, font=2, cex=1.5)
  mtext("神经元凋亡执行通路（caspase/Bim/Bmf）时序激活", side=3, line=0.5, cex=0.95, col=ink)
  legend("topright", legend="* p<0.05 vs NoTBI\nKruskal-Wallis p=0.008",
         bty="n", cex=0.78, text.col=ink)

  ## Panel D: anoikis 倾向 时序
  m  <- tapply(agg$anoikis, agg$timepoint, mean)
  se <- tapply(agg$anoikis, agg$timepoint, function(x) sd(x)/sqrt(length(x)))
  bp <- barplot(m, col=c(navy,steel,accent,accent), border=NA,
                ylim=c(min(m-se)*1.3, max(m+se)+0.06),
                ylab="失巢凋亡倾向指数  (death − survival)", main="")
  abline(h=0, col="#9aa7b3", lty=2)
  arrows(bp, m-se, bp, m+se, angle=90, code=3, length=0.05, col=ink, lwd=1.4)
  set.seed(1)
  for (i in seq_along(levels(agg$timepoint))) {
    tp <- levels(agg$timepoint)[i]; y <- agg$anoikis[agg$timepoint==tp]
    points(jitter(rep(bp[i],length(y)),amount=0.10), y, pch=21, bg="white", col=ink, cex=0.95)
  }
  text(bp[3], m[3]+se[3]+0.020, "*", cex=1.9, col=accent, font=2)
  text(bp[4], m[4]+se[4]+0.020, "*", cex=1.9, col=accent, font=2)
  mtext("D", side=3, line=1.6, adj=-0.20, font=2, cex=1.5)
  mtext("神经元失巢凋亡倾向指数（慢性期显著升高）", side=3, line=0.5, cex=0.95, col=ink)
  legend("topright", legend="* p<0.05 vs NoTBI\nKruskal-Wallis p=0.018",
         bty="n", cex=0.78, text.col=ink)
}

png(file.path(OUT,"fig7_panel_CD.png"), width=2400, height=1050, res=300,
    type="cairo", family="Microsoft YaHei")
draw(); dev.off()
cat("图7 panel C/D 已生成 ->", OUT, "\n")
