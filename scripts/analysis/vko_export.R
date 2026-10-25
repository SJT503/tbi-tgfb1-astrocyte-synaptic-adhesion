#!/usr/bin/env Rscript
# 虚拟敲除 step1 —— 导出整合小鼠星形胶质数据(可靠数据,纯度11.1)
suppressMessages({library(Seurat); library(Matrix)})
P <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project"
OUT <- file.path(P,"JNI/analysis/vko"); dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

obj <- readRDS(file.path(P,"processed_data/integrated/mouse_integrated_annotated.rds"))
obj <- JoinLayers(obj, assay="RNA")
ast <- obj[, obj$celltype_auto=="astrocyte"]
cat("星形胶质细胞:", ncol(ast), "\n")

dat <- GetAssayData(ast, assay="RNA", layer="data")    # 已 log-normalized
# 纯度门禁(铁律)
gi <- function(g){ h<-rownames(dat)[toupper(rownames(dat))==toupper(g)]; if(length(h)) h[1] else NA }
pr <- mean(dat[gi("Aqp4"),])/max(mean(dat[gi("Snap25"),]),1e-6)
cat(sprintf("纯度 Aqp4/Snap25 = %.1f\n", pr))
if (pr < 3) stop("PURITY FAIL")

# 只保留有表达的基因(在>=1%细胞检出),压缩规模
det <- Matrix::rowMeans(dat>0)
keep <- det >= 0.01
dat <- dat[keep,]
cat("保留基因(检出>=1%):", nrow(dat), "\n")

writeMM(dat, file.path(OUT,"data.mtx"))
writeLines(rownames(dat), file.path(OUT,"genes.txt"))
writeLines(colnames(dat), file.path(OUT,"cells.txt"))
md <- data.frame(cell=colnames(ast),
                 condition=as.character(ast$condition),
                 region=as.character(ast$brain_region))
write.csv(md, file.path(OUT,"meta.csv"), row.names=FALSE)
cat("DONE -> ", OUT, "\n")
cat("分组:\n"); print(table(md$condition))
