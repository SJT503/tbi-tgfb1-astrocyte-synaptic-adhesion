#!/usr/bin/env Rscript
# GSE180862 按脑区拆分 (Hip vs Cortex), 取海马星形胶质算 FOXG1 lfc + rho + VKO 位移
suppressMessages({library(Seurat); library(Matrix)})
set.seed(1)
P<-"D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project"
o<-readRDS(file.path(P,"processed_data/integrated/mouse_integrated_annotated.rds"))
o<-JoinLayers(o, assay="RNA")
ast<-o[, o$celltype_auto=="astrocyte" & o$dataset=="GSE180862"]
cat("GSE180862 astrocyte 总:", ncol(ast), "  brain_region:\n"); print(table(ast$brain_region))

ADH<-c("Nrxn1","Nrxn3","Cadm1","Cadm2","Cntn1","Cdh2")
gci<-function(rn,g){h<-rn[toupper(rn)==toupper(g)]; if(length(h)) h[1] else NA}

run_subset<-function(tag, sub){
  d<-GetAssayData(sub,assay="RNA",layer="data")
  rn<-rownames(d)
  pr<-mean(d[gci(rn,"Aqp4"),])/max(mean(d[gci(rn,"Snap25"),]),1e-6)
  cat(sprintf("\n## %s  n=%d  纯度=%.1f\n",tag,ncol(sub),pr))
  foxg<-gci(rn,"Foxg1")
  cond<-as.character(sub$condition); inj<-cond=="tbi"
  cat(sprintf("  损伤态 %d  对照 %d\n", sum(inj), sum(!inj)))
  adg<-na.omit(sapply(ADH,function(g) gci(rn,g)))
  adh_mod<-Matrix::colMeans(d[adg,,drop=FALSE])
  lfc<-log2((mean(d[foxg,inj])+0.01)/(mean(d[foxg,!inj])+0.01))
  rho<-suppressWarnings(cor(d[foxg,inj], adh_mod[inj], method="spearman"))
  cat(sprintf("  FOXG1: lfc=%+.3f  rho_adh=%+.3f\n", lfc, rho))
  # VKO
  det<-Matrix::rowMeans(d>0)>=0.05; gu<-rn[det]
  cells<-1:ncol(sub); if(length(cells)>8000) cells<-sample(cells,8000)
  X<-t(as.matrix(d[gu,cells,drop=FALSE])); mu<-colMeans(X); sd<-apply(X,2,sd); sd[sd==0]<-1
  Z<-scale(X,center=mu,scale=sd)
  tf_mm<-toupper(readLines(file.path(P,"scenic_input/databases/allTFs_mm.txt")))
  tfs<-intersect(gu, rn[toupper(rn)%in%tf_mm])
  tv<-apply(Z[,tfs,drop=FALSE],2,var); tfs<-names(sort(tv,decreasing=TRUE))[1:min(150,length(tfs))]
  if(!(foxg %in% tfs)) tfs<-c(tfs[-length(tfs)],foxg)
  Xtf<-Z[,tfs,drop=FALSE]
  B<-solve(crossprod(Xtf)+diag(length(tfs)), crossprod(Xtf,Z))
  fp<-which(tfs==foxg); dF<-((0-mu[foxg])/sd[foxg])-Z[,foxg]
  s1<-outer(dF,B[fp,]); s2<-s1+s1[,tfs,drop=FALSE]%*%B
  ms<-colMeans(s2); names(ms)<-colnames(Z)
  ash<-ms[adg]
  u<-suppressWarnings(wilcox.test(ash,ms,alternative="less"))
  cat(sprintf("  VKO 黏附均值=%+.4f  下调 %d/%d  p=%.3f\n",
      mean(ash,na.rm=T),sum(ash<0,na.rm=T),length(adg),u$p.value))
  invisible(list(tag=tag,n=ncol(sub),pr=pr,lfc=round(lfc,3),rho=round(rho,3),
                 ks=round(mean(ash,na.rm=T),4),
                 down=sprintf("%d/%d",sum(ash<0,na.rm=T),length(adg)),
                 p=signif(u$p.value,3)))
}

r_hip<-run_subset("GSE180862-Hip", ast[, ast$brain_region=="hippocampus"])
r_cor<-run_subset("GSE180862-Cortex", ast[, ast$brain_region=="cortex"])

cat("\n\n=== 汇总 ===\n")
cat(sprintf("GSE180862-Hip    n=%d pr=%.1f  lfc=%+.3f rho=%+.3f  VKO=%+.4f %s p=%.3f\n",
  r_hip$n,r_hip$pr,r_hip$lfc,r_hip$rho,r_hip$ks,r_hip$down,r_hip$p))
cat(sprintf("GSE180862-Cortex n=%d pr=%.1f  lfc=%+.3f rho=%+.3f  VKO=%+.4f %s p=%.3f\n",
  r_cor$n,r_cor$pr,r_cor$lfc,r_cor$rho,r_cor$ks,r_cor$down,r_cor$p))
