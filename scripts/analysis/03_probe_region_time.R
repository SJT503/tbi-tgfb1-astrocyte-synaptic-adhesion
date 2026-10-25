#!/usr/bin/env Rscript
suppressMessages({library(Seurat)})
seu <- readRDS("D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/processed_data/integrated/mouse_integrated_annotated.rds")
mt <- seu@meta.data

cat("===== dataset x brain_region =====\n")
print(table(mt$dataset, mt$brain_region))

cat("\n===== dataset x brain_region x timepoint =====\n")
print(xtabs(~ dataset + brain_region + timepoint, data=mt))

cat("\n===== dataset x brain_region x timepoint x condition =====\n")
print(xtabs(~ dataset + brain_region + timepoint + condition, data=mt))

cat("\n===== animal x brain_region x timepoint x condition =====\n")
print(xtabs(~ animal + brain_region + timepoint + condition, data=mt))

cat("\n===== celltype_auto x brain_region =====\n")
print(table(mt$celltype_auto, mt$brain_region))

cat("\n===== celltype_auto x brain_region x condition (NEURON ONLY) =====\n")
neu_mt <- mt[mt$celltype_auto == "neuron", ]
print(xtabs(~ brain_region + timepoint + condition, data=neu_mt))

cat("\n===== mice biology replication (per brain region x timepoint x condition) =====\n")
neu_mt$key <- paste(neu_mt$brain_region, neu_mt$timepoint, neu_mt$condition, sep="|")
agg <- aggregate(animal ~ key, data=neu_mt, FUN=function(x) length(unique(x)))
print(agg)
