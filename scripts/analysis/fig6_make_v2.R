#!/usr/bin/env Rscript
# 研究基础图6 v2 — 按批注#4修订:
# B:用 GSE180862-Hip 替换效应弱的 GSE101901; 顺序调换使人源 GSE209552 居右并高亮
# C:同上, GSE226207 因纯度不合格剔除
# D:重排 60%/69% 标签与 legend 位置, 消除重叠
# X轴 GSE226207 标签从"纯度1.6*" 改为"(纯度不合格)" 避免歧义
A<-"D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis"
OUT<-file.path(A,"fig6_out"); dir.create(OUT,showWarnings=FALSE)
navy<-"#1a3a5c";mid<-"#3f6896";steel<-"#7795b5";accent<-"#b5563f";teal<-"#2a6a5e"
grey1<-"#c8d2dc";ink<-"#2a3a4a";gold<-"#c8a45c"

scr<-read.csv(file.path(A,"tf_screen_out/GSE101901_180862_ms_TFscreen.csv"))
tc<-read.csv(file.path(A,"foxg1_timecourse/foxg1_astro_per_sample.csv"))

# 5 数据集新顺序(人源在右,高亮)
v<-data.frame(
 short=c("GSE180862-Hip","GSE180862-Ctx","GSE226207","GSE271769","GSE209552"),
 lab  =c("GSE180862\n海马","GSE180862\n皮层","GSE226207\n(纯度不合格)","GSE271769\n重度CCI","GSE209552\n人源·Human"),
 lfc  =c(-0.285,-0.234,-0.237,-0.510,-0.760),
 rho  =c(+0.041,+0.081,+0.275,+0.205,+0.244),
 ko   =c(-0.0025,-0.0024,NA,-0.0198,-0.0163),
 p    =c(0.392,0.366,NA,0.027,0.074),
 hum  =c(FALSE,FALSE,FALSE,FALSE,TRUE),
 excl =c(FALSE,FALSE,TRUE,FALSE,FALSE)
)
col_b<-ifelse(v$excl,grey1,ifelse(v$hum,accent,mid))
col_c<-ifelse(v$excl,grey1,ifelse(v$hum,accent,ifelse(v$p<0.05,navy,steel)))

draw<-function(){
  par(mfrow=c(2,2), mar=c(5.0,4.8,3.4,1.4), mgp=c(2.7,0.7,0), family="Microsoft YaHei",
      cex.axis=0.92, cex.lab=1.05)

  ## A: TF筛选散点
  plot(scr$lfc_injury, scr$rho_within, pch=19, cex=0.7, col=grey1,
       xlab="损伤效应  log2FC（损伤 vs 对照）",
       ylab="与突触黏附程序相关  Spearman r", main="")
  abline(h=0,v=0,lty=2,col="#9aa7b3")
  fx<-scr[toupper(scr$TF)=="FOXG1",][1,]
  points(fx$lfc_injury, fx$rho_within, pch=21, bg=accent, col="white", cex=2.3, lwd=2)
  text(fx$lfc_injury, fx$rho_within, "FOXG1", pos=2, font=2, col=accent, cex=1.05, offset=0.7)
  mtext("A",side=3,line=1.6,adj=-0.18,font=2,cex=1.5)
  mtext("207 个转录因子的损伤效应×黏附相关分布",side=3,line=0.5,cex=0.92,col=ink)
  legend("topright", legend="损伤下调 + 与黏附正相关\n= 丢失的黏附正调控者",
         bty="n", text.col=ink, cex=0.78)

  ## B: FOXG1 损伤后下调, 5 个独立数据集(人源最右且高亮)
  bp<-barplot(v$lfc, names.arg=NA, col=col_b, border=NA,
       ylab="FOXG1 损伤效应  log2FC",
       ylim=c(min(v$lfc)*1.30,0.18), main="")
  abline(h=0,col="#9aa7b3")
  text(bp, v$lfc-0.045, sprintf("%.2f",v$lfc), cex=0.84, col=ink, font=2)
  text(bp, rep(0.10,5), sprintf("r=%.2f",v$rho), cex=0.78, col=teal)
  text(bp, par("usr")[3]-0.05, v$lab, xpd=TRUE, cex=0.74, col=ink)
  # 人源加边框强调
  hidx<-which(v$hum)
  rect(bp[hidx]-0.5, par("usr")[3], bp[hidx]+0.5, par("usr")[4],
       border=accent, lwd=2.2, lty=3)
  mtext("B",side=3,line=1.6,adj=-0.18,font=2,cex=1.5)
  mtext("FOXG1 在 5 个独立数据集损伤后一致下调（人源居右）",side=3,line=0.5,cex=0.92,col=ink)

  ## C: 虚拟敲除黏附程序位移
  vc<-v[!is.na(v$ko),]
  col_cc<-ifelse(vc$excl,grey1,ifelse(vc$hum,accent,ifelse(vc$p<0.05,navy,steel)))
  bp2<-barplot(vc$ko*1000, names.arg=NA, col=col_cc, border=NA,
       ylab="虚拟敲除黏附程序位移  (×0.001)",
       ylim=c(min(vc$ko*1000,na.rm=T)*1.35, 4), main="")
  abline(h=0,col="#9aa7b3")
  sig<-ifelse(vc$p<0.05,"*",ifelse(vc$p<0.1,"#",""))
  text(bp2, vc$ko*1000-1.5, sprintf("p=%.3f",vc$p), cex=0.78, col=ink)
  text(bp2, vc$ko*1000-3.0, sig, cex=1.7, col=accent, font=2)
  c_lab<-c("GSE180862\nHip","GSE180862\nCtx","GSE271769\n重度CCI","GSE209552\n人源")
  text(bp2, par("usr")[3]-1.5, c_lab[match(vc$short,c("GSE180862-Hip","GSE180862-Ctx","GSE271769","GSE209552"))],
       xpd=TRUE, cex=0.72, col=ink)
  hidx2<-which(vc$hum)
  rect(bp2[hidx2]-0.5, par("usr")[3], bp2[hidx2]+0.5, par("usr")[4],
       border=accent, lwd=2.2, lty=3)
  mtext("C",side=3,line=1.6,adj=-0.18,font=2,cex=1.5)
  mtext("FOXG1 虚拟敲除预测突触黏附程序下移",side=3,line=0.5,cex=0.92,col=ink)
  legend("bottomleft", legend=c("p<0.05","p<0.1","人源数据"),
         fill=c(navy,steel,accent), border=NA, bty="n", cex=0.74)

  ## D: FOXG1 时序(GSE271769), 修正标签重叠
  tc$timepoint<-factor(tc$timepoint, levels=c("NoTBI","D1","D7","D28"))
  m <-tapply(tc$foxg1, tc$timepoint, mean)
  se<-tapply(tc$foxg1, tc$timepoint, function(x) sd(x)/sqrt(length(x)))
  bp3<-barplot(m, col=c(navy,accent,steel,accent), border=NA, ylim=c(0,1.15),
       ylab="FOXG1 星形胶质细胞表达量", main="")
  arrows(bp3, m-se, bp3, m+se, angle=90, code=3, length=0.05, col=ink, lwd=1.4)
  set.seed(1)
  for(i in seq_along(levels(tc$timepoint))){
    tp<-levels(tc$timepoint)[i]; y<-tc$foxg1[tc$timepoint==tp]
    points(jitter(rep(bp3[i],length(y)),amount=0.09), y, pch=21, bg="white", col=ink, cex=0.9)
  }
  # 显著性*与百分比标签紧贴误差棒上方,不再放在固定高位
  pct<-c("","60%","","69%")
  star<-c("","*","","*")
  for(i in 1:4){
    if(star[i]!="" ){
      text(bp3[i], m[i]+se[i]+0.06, star[i], cex=1.7, col=accent, font=2)
      text(bp3[i], m[i]+se[i]+0.13, pct[i], cex=0.82, col=ink, font=2)
    }
  }
  mtext("D",side=3,line=1.6,adj=-0.18,font=2,cex=1.5)
  mtext("FOXG1 在 TBI 后皮层星形胶质细胞的时序表达（GSE271769）",side=3,line=0.5,cex=0.88,col=ink)
  # legend 移到底部中间, 避开柱子顶端
  legend("bottom", inset=c(0,-0.30), horiz=TRUE,
         legend=c("* p<0.05 vs NoTBI","Kruskal-Wallis p=0.021"),
         bty="n", cex=0.78, text.col=ink, xpd=TRUE)
}

png(file.path(OUT,"研究基础图6_v2.png"), width=2600, height=2180, res=300, type="cairo", family="Microsoft YaHei")
draw(); dev.off()
cairo_pdf(file.path(OUT,"研究基础图6_v2.pdf"), width=8.7, height=7.3, family="Microsoft YaHei"); draw(); dev.off()
cat("图6 v2 已生成 ->",OUT,"\n")
