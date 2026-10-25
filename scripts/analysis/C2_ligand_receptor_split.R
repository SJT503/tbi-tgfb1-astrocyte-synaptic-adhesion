#!/usr/bin/env Rscript
# C2: Ligand (astrocyte side) vs Receptor (neuron side) split — 精准锁定塌陷驱动端
# 从已跑的 CellChat RDS 提:
#   (a) 公共小鼠 cand2: control vs tbi
#   (b) 人源 cand2b: control vs tbi
#   (c) 自测三组 wet_lab: Sham/CCI_Scr/CCI_KO
# 输出: ligand 端 expression 变化 vs receptor 端 expression 变化 (NRXN family)
suppressMessages({library(CellChat); library(ggplot2); library(patchwork); library(dplyr)})

OUT <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/C2_ligand_receptor"
WET <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/wet_lab_data/results/cellchat"
PUB <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/cand_out"

# adhesion ligands (astrocyte side) and receptors (neuron side)
LIGANDS  <- c("Nrxn1","Nrxn2","Nrxn3","Cadm1","Cadm2","Cadm3","Cntn1","Cntn2","Cntn4")
RECEPTORS<- c("Nlgn1","Nlgn2","Nlgn3","Cadm1","Cadm2","Cadm3","Cntn1","Cntn2",
              "Lrrtm1","Lrrtm2","Lrrtm3","Lrrtm4")

extract_cell_expr <- function(cc, gene_list, celltype_pattern) {
  # cc: CellChat object; gene_list: genes to extract; celltype_pattern: regex
  data_use <- cc@data.signaling  # genes × cells, normalized
  meta <- cc@meta
  ct_col <- if ("labels" %in% colnames(meta)) "labels" else colnames(meta)[1]
  ct_vec <- as.character(meta[[ct_col]])
  cells_hit <- grepl(celltype_pattern, ct_vec, ignore.case=TRUE)
  if (sum(cells_hit) == 0) return(NULL)
  res <- sapply(gene_list, function(g) {
    if (!g %in% rownames(data_use)) return(NA)
    mean(data_use[g, cells_hit])
  })
  data.frame(gene=names(res), mean_expr=unname(res),
             n_cells=sum(cells_hit), celltype_used=ct_col)
}

cat("[C2.1] PUBLIC mouse cand2 (control vs tbi)\n")
res_pub <- list()
for (cond in c("control","tbi")) {
  f <- file.path(PUB, sprintf("cand2_cellchat_%s.rds", cond))
  if (!file.exists(f)) { cat("  missing:", f, "\n"); next }
  cc <- readRDS(f)
  cat("  loaded", cond, ": cells=", ncol(cc@data.signaling),
      " celltypes=", paste(levels(cc@meta$labels), collapse=","), "\n")
  lig_astro <- extract_cell_expr(cc, LIGANDS,   "astro|Ast")
  rec_neuro <- extract_cell_expr(cc, RECEPTORS, "neuron|Neur|Excit|Inhib")
  if (!is.null(lig_astro)) { lig_astro$side <- "ligand_on_astrocyte"; lig_astro$group <- cond }
  if (!is.null(rec_neuro)) { rec_neuro$side <- "receptor_on_neuron"; rec_neuro$group <- cond }
  res_pub[[cond]] <- rbind(lig_astro, rec_neuro)
}
pub_df <- do.call(rbind, res_pub)
if (!is.null(pub_df)) {
  write.csv(pub_df, file.path(OUT, "C2_public_mouse_LR.csv"), row.names=FALSE)
  cat("  saved C2_public_mouse_LR.csv\n")
}

cat("\n[C2.2] PUBLIC human cand2b\n")
res_hum <- list()
for (cond in c("control","tbi")) {
  f <- file.path(PUB, sprintf("cand2b_cellchat_%s.rds", cond))
  if (!file.exists(f)) next
  cc <- readRDS(f)
  L_h <- toupper(LIGANDS); R_h <- toupper(RECEPTORS)
  lig <- extract_cell_expr(cc, L_h, "astro|Ast")
  rec <- extract_cell_expr(cc, R_h, "neuron|Neur|Excit|Inhib")
  if (!is.null(lig)) { lig$side <- "ligand_on_astrocyte"; lig$group <- cond }
  if (!is.null(rec)) { rec$side <- "receptor_on_neuron"; rec$group <- cond }
  res_hum[[cond]] <- rbind(lig, rec)
}
hum_df <- do.call(rbind, res_hum)
if (!is.null(hum_df)) {
  write.csv(hum_df, file.path(OUT, "C2_human_LR.csv"), row.names=FALSE)
  cat("  saved C2_human_LR.csv\n")
}

cat("\n[C2.3] SELF three-group wet_lab\n")
res_self <- list()
for (grp in c("Sham","CCI_Scr","CCI_TGFbetaR1_KO")) {
  f <- file.path(WET, sprintf("cellchat_%s.rds", grp))
  if (!file.exists(f)) { cat("  missing:", f, "\n"); next }
  cc <- readRDS(f)
  cat("  loaded", grp, ": celltypes=", paste(levels(cc@meta$labels), collapse=","), "\n")
  lig <- extract_cell_expr(cc, LIGANDS,   "astro|Ast")
  rec <- extract_cell_expr(cc, RECEPTORS, "neuron|Neur|DG|CA")
  if (!is.null(lig)) { lig$side <- "ligand_on_astrocyte"; lig$group <- grp }
  if (!is.null(rec)) { rec$side <- "receptor_on_neuron"; rec$group <- grp }
  res_self[[grp]] <- rbind(lig, rec)
}
self_df <- do.call(rbind, res_self)
if (!is.null(self_df)) {
  write.csv(self_df, file.path(OUT, "C2_self_LR.csv"), row.names=FALSE)
  cat("  saved C2_self_LR.csv\n")
}

cat("\n[C2.4] driver-side summary\n")
# For each dataset, which side drops more after TBI?
summary_driver <- function(df, label) {
  if (is.null(df)) return(NULL)
  df %>% tidyr::pivot_wider(id_cols=c(gene, side), names_from=group, values_from=mean_expr) -> w
  cat("\n", label, "\n", sep="")
  print(w)
  w
}
summary_driver(pub_df,  "PUBLIC mouse cand2")
summary_driver(hum_df,  "PUBLIC human cand2b")
summary_driver(self_df, "SELF three-group")

# ── FIGURES ─────────────────────────────────────────────────────────
make_panel <- function(df, title) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  ggplot(df, aes(gene, mean_expr, fill=group)) +
    geom_col(position=position_dodge(0.8), width=0.7, colour="grey25", linewidth=0.2) +
    facet_wrap(~side, scales="free", ncol=1) +
    labs(title=title, y="mean normalized expression", x=NULL) +
    theme_bw(base_size=10) +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top")
}
p_pub  <- make_panel(pub_df,  "C2A. Public mouse — ligand (astro) vs receptor (neuron)")
p_hum  <- make_panel(hum_df,  "C2B. Human — ligand (astro) vs receptor (neuron)")
p_self <- make_panel(self_df, "C2C. Self three-group — ligand vs receptor (descriptive)")

for (pp in list(p_pub, p_hum, p_self)) {
  if (is.null(pp)) next
  ttl <- gsub("[.: —]", "_", pp$labels$title)
  ggsave(file.path(OUT, paste0(substr(ttl, 1, 20), ".pdf")), pp, width=8, height=5)
  ggsave(file.path(OUT, paste0(substr(ttl, 1, 20), ".png")), pp, width=8, height=5, dpi=200)
}

cat("\n[C2.5] DONE at:\n  ", OUT, "\n")
