#!/usr/bin/env Rscript
# 研究基础图6（发表级，4 面板）—— TF 筛选锁定 FOXG1 + 跨数据集验证 + 虚拟敲除 + 时序
# 数据：tf_screen_out / vko_per_dataset / foxg1_timecourse  (均为真实分析输出)
A <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis"
OUT <- file.path(A,"fig6_out"); dir.create(OUT, showWarnings=FALSE)

navy<-"#1a3a5c"; mid<-"#3f6896"; steel<-"#7795b5"; accent<-"#b5563f"; teal<-"#2a6a5e"
grey1<-"#c8d2dc"; ink<-"#2a3a4a"

scr <- read.csv(file.path(A,"tf_screen_out/GSE101901_180862_ms_TFscreen.csv"))
vko <- read.csv(file.path(A,"vko_per_dataset/vko_per_dataset_summary.csv"))
tc  <- read.csv(file.path(A,"foxg1_timecourse/foxg1_astro_per_sample.csv"))

draw <- function(){
  par(mfrow=c(2,2), mar=c(4.6,4.8,3.4,1.4), mgp=c(2.7,0.7,0), family="sans",
      cex.axis=0.92, cex.lab=1.05)

  ## ---- A：TF 筛选全景散点 ----
  plot(scr$lfc_injury, scr$rho_within, pch=19, cex=0.7, col=grey1,
       xlab="损伤效应  log2FC（损伤 vs 对照）", ylab="与突触黏附程序相关  Spearman ρ",
       main="")
  abline(h=0,v=0,lty=2,col="#9aa7b3")
  fx <- scr[toupper(scr$TF)=="FOXG1",][1,]
  points(fx$lfc_injury, fx$rho_within, pch=21, bg=accent, col="white", cex=2.3, lwd=2)
  text(fx$lfc_injury, fx$rho_within, "FOXG1", pos=2, font=2, col=accent, cex=1.05, offset=0.7)
  mtext("A", side=3, line=1.6, adj=-0.18, font=2, cex=1.5)
  mtext("207 个转录因子的损伤效应×黏附相关分布", side=3, line=0.5, cex=0.92, col=ink)
  legend("topright", legend="损伤下调 + 与黏附正相关\n= 丢失的黏附正调控者", bty="n",
         text.col=ink, cex=0.78)

  ## ---- B：FOXG1 损伤后下调，5 个独立数据集（拆开） ----
  dord <- c("GSE101901","GSE180862","GSE226207","GSE209552","GSE271769")
  v <- vko[match(dord, vko$dataset),]
  cols <- ifelse(v$dataset=="GSE226207", grey1, mid)
  bp <- barplot(v$foxg1_lfc_injury, names.arg=NA, col=cols, border=NA,
       ylab="FOXG1 损伤效应  log2FC", ylim=c(min(v$foxg1_lfc_injury)*1.25,0.12), main="")
  abline(h=0,col="#9aa7b3")
  text(bp, v$foxg1_lfc_injury-0.04, sprintf("%.2f",v$foxg1_lfc_injury), cex=0.82, col=ink)
  text(bp, rep(0.05,5), sprintf("ρ=%.2f",v$foxg1_rho_adh), cex=0.74, col=teal)
  text(bp, par("usr")[3]-0.04, c("GSE101901\n海马","GSE180862\nFPI","GSE226207\n纯度1.6*","GSE209552\n人源","GSE271769\n重度CCI"),
       xpd=TRUE, cex=0.72, col=ink, srt=0)
  mtext("B", side=3, line=1.6, adj=-0.18, font=2, cex=1.5)
  mtext("FOXG1 在 5 个独立数据集损伤后一致下调", side=3, line=0.5, cex=0.92, col=ink)

  ## ---- C：FOXG1 虚拟敲除 → 黏附程序预测位移（跨数据集） ----
  vc <- v[!is.na(v$ko_adh_mean_shift),]
  cols2 <- ifelse(vc$dataset=="GSE226207", grey1,
            ifelse(vc$ko_wilcox_p<0.05, accent, steel))
  bp2 <- barplot(vc$ko_adh_mean_shift*1000, names.arg=NA, col=cols2, border=NA,
       ylab="虚拟敲除黏附程序位移  (×10⁻³)", ylim=c(min(vc$ko_adh_mean_shift*1000)*1.3,3), main="")
  abline(h=0,col="#9aa7b3")
  sig <- ifelse(vc$ko_wilcox_p<0.05,"*", ifelse(vc$ko_wilcox_p<0.1,"#",""))
  text(bp2, vc$ko_adh_mean_shift*1000-1.2, sprintf("p=%.3f",vc$ko_wilcox_p), cex=0.74, col=ink)
  text(bp2, vc$ko_adh_mean_shift*1000-2.6, sig, cex=1.5, col=accent)
  text(bp2, par("usr")[3]-1.1, vc$dataset, xpd=TRUE, cex=0.74, col=ink)
  mtext("C", side=3, line=1.6, adj=-0.18, font=2, cex=1.5)
  mtext("FOXG1 虚拟敲除预测突触黏附程序下移", side=3, line=0.5, cex=0.92, col=ink)
  legend("bottomleft", legend=c("p<0.05","p<0.1","纯度不合格(剔除)"),
         fill=c(accent,steel,grey1), border=NA, bty="n", cex=0.72)

  ## ---- D：FOXG1 时序表达（GSE271769） ----
  tc$timepoint <- factor(tc$timepoint, levels=c("NoTBI","D1","D7","D28"))
  m  <- tapply(tc$foxg1, tc$timepoint, mean)
  se <- tapply(tc$foxg1, tc$timepoint, function(x) sd(x)/sqrt(length(x)))
  bp3 <- barplot(m, col=c(navy,accent,steel,accent), border=NA, ylim=c(0,1.0),
       ylab="FOXG1 星形胶质细胞表达量", main="")
  arrows(bp3, m-se, bp3, m+se, angle=90, code=3, length=0.05, col=ink, lwd=1.4)
  set.seed(1)
  for(i in seq_along(levels(tc$timepoint))){
    tp<-levels(tc$timepoint)[i]; y<-tc$foxg1[tc$timepoint==tp]
    points(jitter(rep(bp3[i],length(y)),amount=0.09), y, pch=21, bg="white", col=ink, cex=0.9)
  }
  text(bp3[2], m[2]+se[2]+0.07, "*", cex=1.7, col=accent)
  text(bp3[4], m[4]+se[4]+0.07, "*", cex=1.7, col=accent)
  text(bp3[2], 0.93, "60%", cex=0.8, col=ink); text(bp3[4], 0.93, "69%", cex=0.8, col=ink)
  mtext("D", side=3, line=1.6, adj=-0.18, font=2, cex=1.5)
  mtext("FOXG1 在 TBI 后皮层星形胶质细胞的时序表达", side=3, line=0.5, cex=0.92, col=ink)
  legend("topright", legend="* p<0.05 vs NoTBI\nKruskal-Wallis p=0.021", bty="n", cex=0.74, text.col=ink)
}

png(file.path(OUT,"研究基础图6.png"), width=2600, height=2150, res=300)
draw(); dev.off()
pdf(file.path(OUT,"研究基础图6.pdf"), width=8.7, height=7.2)
draw(); dev.off()
cat("图6 已生成 ->", OUT, "\n")
