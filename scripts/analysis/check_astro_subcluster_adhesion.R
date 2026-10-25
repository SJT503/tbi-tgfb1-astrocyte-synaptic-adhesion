#!/usr/bin/env Rscript
# 亚群层次拆解: Astro_C0/C1/C2 → Neuron 在 NRXN/NCAM/NGL 上各组方向
suppressMessages({library(CellChat)})
D <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/wet_lab_data/results/cellchat_nrxn"
cc <- list(
  Sham = readRDS(file.path(D,"cellchat_nrxn_Sham.rds")),
  CCI  = readRDS(file.path(D,"cellchat_nrxn_CCI_Scr.rds")),
  KO   = readRDS(file.path(D,"cellchat_nrxn_CCI_TGFbetaR1_KO.rds"))
)

ASTROS <- c("Astro_C0","Astro_C1","Astro_C2")
NEURONS <- c("Excitatory_Neuron","Inhibitory_Neuron")
PATHS <- c("NRXN","NCAM","NGL","NEGR")

# 4D 矩阵: pathway × astro_sub × neuron × group
results <- list()
for (pn in PATHS){
  cat("\n=== PATHWAY:", pn, "===\n")
  cat(sprintf("%-12s %-22s %10s %10s %10s\n", "AstroSub", "Neuron", "Sham", "CCI", "CCI+KO"))
  for (a in ASTROS){
    for (n in NEURONS){
      vals <- numeric(3)
      for (gi in seq_along(cc)){
        obj <- cc[[gi]]
        idents <- levels(obj@idents)
        if (! pn %in% obj@netP$pathways) { vals[gi] <- NA; next }
        pi <- which(obj@netP$pathways == pn)
        if (! a %in% idents || ! n %in% idents) { vals[gi] <- NA; next }
        ai <- which(idents==a); ni <- which(idents==n)
        vals[gi] <- obj@netP$prob[ai, ni, pi]
      }
      cat(sprintf("%-12s %-22s %10.4f %10.4f %10.4f", a, n, vals[1], vals[2], vals[3]))
      # 显示百分比 vs Sham
      if (!is.na(vals[1]) && vals[1] > 0)
        cat(sprintf("    | %.0f%% / %.0f%% / %.0f%%",
                    100, vals[2]/vals[1]*100, vals[3]/vals[1]*100))
      cat("\n")
      results[[paste(pn,a,n,sep="|")]] <- vals
    }
  }
}
saveRDS(results, "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/astro_sub_adhesion_3groups.rds")
cat("\n done\n")
