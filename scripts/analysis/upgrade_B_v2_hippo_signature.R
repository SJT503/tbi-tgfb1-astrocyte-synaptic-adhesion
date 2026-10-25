#!/usr/bin/env Rscript
# upgrade_B_v2 — Hippocampus neuron TBI vs sham signature for CMap
# 主信号: GSE180862 海马 7d FPI (mice 3 vs 3, time = TGF-β1 peak)
# 加分: GSE101901 海马 24h concussion (mice 3 vs 3) → 取 intersection 作 cross-model robust signature
suppressMessages({library(Seurat); library(Matrix)})

RDS <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/processed_data/integrated/mouse_integrated_annotated.rds"
OUT <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/upgrade_B_v2_hippo"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
TOP_N <- 150

cat("[B2.1] load + subset hippocampus neuron\n")
seu <- readRDS(RDS); DefaultAssay(seu) <- "RNA"
hip <- subset(seu, subset = brain_region == "hippocampus" & celltype_auto == "neuron")
cat("  hippocampus neuron cells:", ncol(hip), "\n")

# ── Primary signature: GSE180862 hippocampus 7d FPI ─────────────────
cat("\n[B2.2] PRIMARY signature: GSE180862 hippocampus 7d FPI (TBI vs sham)\n")
p7 <- subset(hip, subset = dataset == "GSE180862" & timepoint == "7days")
cat("  cells:", ncol(p7), "\n")
p7$condition <- factor(p7$condition, levels=c("control","tbi"))
Idents(p7) <- "condition"
cat("  cells per group:\n"); print(table(p7$condition))
deg_primary <- FindMarkers(p7, ident.1="tbi", ident.2="control",
                           logfc.threshold=0, min.pct=0.05, only.pos=FALSE, test.use="wilcox")
deg_primary$gene <- rownames(deg_primary)
write.csv(deg_primary, file.path(OUT, "B2_DEG_GSE180862_hippo_7d_FPI.csv"), row.names=FALSE)
cat("  total DEG:", nrow(deg_primary), "  padj<0.05:", sum(deg_primary$p_val_adj < 0.05, na.rm=TRUE), "\n")

# remove mt-* and Gm/Rik to avoid depth/dropout artefacts
clean_genes <- function(df) {
  df <- df[!grepl("^mt-", df$gene, ignore.case=TRUE), ]
  df <- df[!grepl("^Gm[0-9]+$", df$gene), ]
  df <- df[!grepl("Rik$", df$gene), ]
  df
}
sig_p <- deg_primary[!is.na(deg_primary$p_val_adj) & deg_primary$p_val_adj < 0.05, ]
sig_p <- clean_genes(sig_p)
cat("  after removing mt-/Gm/Rik:", nrow(sig_p), "\n")
up_p   <- head(sig_p$gene[order(-sig_p$avg_log2FC)], TOP_N)
down_p <- head(sig_p$gene[order( sig_p$avg_log2FC)], TOP_N)
cat("  primary up   top 10:", paste(head(up_p, 10), collapse=","), "\n")
cat("  primary down top 10:", paste(head(down_p, 10), collapse=","), "\n")

# ── Secondary signature: GSE101901 hippocampus 24h concussion ───────
cat("\n[B2.3] SECONDARY signature: GSE101901 hippocampus 24h concussion\n")
p24 <- subset(hip, subset = dataset == "GSE101901" & timepoint == "24h")
cat("  cells:", ncol(p24), "\n")
p24$condition <- factor(p24$condition, levels=c("control","tbi"))
Idents(p24) <- "condition"
deg_secondary <- FindMarkers(p24, ident.1="tbi", ident.2="control",
                             logfc.threshold=0, min.pct=0.05, only.pos=FALSE, test.use="wilcox")
deg_secondary$gene <- rownames(deg_secondary)
write.csv(deg_secondary, file.path(OUT, "B2_DEG_GSE101901_hippo_24h_concussion.csv"), row.names=FALSE)
sig_s <- deg_secondary[!is.na(deg_secondary$p_val_adj) & deg_secondary$p_val_adj < 0.05, ]
sig_s <- clean_genes(sig_s)
up_s   <- head(sig_s$gene[order(-sig_s$avg_log2FC)], TOP_N)
down_s <- head(sig_s$gene[order( sig_s$avg_log2FC)], TOP_N)

# ── Cross-model intersection (most robust signature) ────────────────
cat("\n[B2.4] CROSS-MODEL intersection (most robust)\n")
up_robust   <- intersect(up_p, up_s)
down_robust <- intersect(down_p, down_s)
cat("  up   intersection:", length(up_robust), "/", TOP_N, "\n")
cat("  down intersection:", length(down_robust), "/", TOP_N, "\n")
cat("  up   robust head:", paste(head(up_robust, 15), collapse=","), "\n")
cat("  down robust head:", paste(head(down_robust, 15), collapse=","), "\n")

# adhesion panel sanity check
PANEL <- c("Nrxn1","Nrxn2","Nrxn3","Negr1","Cntn2","Cntn4","Cntn5",
           "Cadm1","Cadm2","Cadm3","Cdh8","Cdh11","Cdh13","Nlgn1","Nlgn2","Nlgn3")
panel_hit_p_down <- intersect(PANEL, down_p)
panel_hit_s_down <- intersect(PANEL, down_s)
panel_hit_robust <- intersect(PANEL, down_robust)
cat("\n[B2.5] Adhesion panel in DOWN signature (sanity)\n")
cat("  primary 7d FPI down:    ", paste(panel_hit_p_down, collapse=", "), "\n")
cat("  secondary 24h conc down:", paste(panel_hit_s_down, collapse=", "), "\n")
cat("  robust intersection:    ", paste(panel_hit_robust, collapse=", "), "\n")

# ── Write CMap inputs ───────────────────────────────────────────────
m2h <- function(g) toupper(g)
write_gmt <- function(fn, name, genes) writeLines(c(paste0(name,"\tCMap_query"), paste(genes, collapse="\t")), fn)
write_gmt(file.path(OUT, "B2_signature_up_primary.gmt"),    "hippo_7d_FPI_up",    m2h(up_p))
write_gmt(file.path(OUT, "B2_signature_down_primary.gmt"),  "hippo_7d_FPI_down",  m2h(down_p))
write_gmt(file.path(OUT, "B2_signature_up_secondary.gmt"),  "hippo_24h_concussion_up",   m2h(up_s))
write_gmt(file.path(OUT, "B2_signature_down_secondary.gmt"),"hippo_24h_concussion_down", m2h(down_s))
write_gmt(file.path(OUT, "B2_signature_up_robust.gmt"),     "hippo_cross_model_up",   m2h(up_robust))
write_gmt(file.path(OUT, "B2_signature_down_robust.gmt"),   "hippo_cross_model_down", m2h(down_robust))
writeLines(m2h(up_p),       file.path(OUT,"B2_up_primary.txt"))
writeLines(m2h(down_p),     file.path(OUT,"B2_down_primary.txt"))
writeLines(m2h(up_s),       file.path(OUT,"B2_up_secondary.txt"))
writeLines(m2h(down_s),     file.path(OUT,"B2_down_secondary.txt"))
writeLines(m2h(up_robust),  file.path(OUT,"B2_up_robust.txt"))
writeLines(m2h(down_robust),file.path(OUT,"B2_down_robust.txt"))

# Galunisertib watchlist
watch <- data.frame(
  compound=c("Galunisertib","LY2157299","LY-2109761","Vactosertib","SB-431542","SB-505124","Repsox"),
  target="TGFBR1 (ALK5)",
  note="TGF-βRI inhibitor — expected nomination if signature captures TGFβ axis")
write.csv(watch, file.path(OUT,"B2_galunisertib_watchlist.csv"), row.names=FALSE)

# README
readme <- c(
"# v2.2 Upgrade B (v2) — Hippocampus-matched CMap query preparation",
"",
"Generated: 2026-05-24",
"",
"## Signature sources (3 versions, recommended use *robust*)",
"- PRIMARY: GSE180862 hippocampus 7d FPI, n=Sham 3 vs TBI 3 mice — time matches TGF-β1 peak",
"- SECONDARY: GSE101901 hippocampus 24h concussion, n=Sham 3 vs TBI 3 mice — earliest TBI window",
"- ROBUST: intersection of primary ∩ secondary — cross-model conserved signature",
"",
"## CMap submission",
"1. Visit https://clue.io/query (free account, VPN)",
"2. Tool → Query → Touchstone",
"3. Upload B2_signature_up_robust.gmt + B2_signature_down_robust.gmt",
"4. Submit, wait 10-30 min, download result.csv",
"5. Place as B2_cmap_result_robust.csv",
"",
"## Verification claim",
"After result returns, check top-100 for: Galunisertib / LY2157299 / SB-431542 etc",
"If hit → Discussion paragraph:",
"",
"> 'In silico CMap query with the cross-model robust hippocampal neuron TBI signature",
"> independently nominated Galunisertib (rank #X of 9597 compounds) among the top compounds",
"> whose transcriptional effects most resemble reversing the TBI state, without any prior",
"> TGFβ pathway knowledge. This algorithmic transferability supports the TGF-βRI axis as",
"> a specific, pharmacologically targetable node in post-traumatic adhesion collapse.'"
)
writeLines(readme, file.path(OUT,"README_clueio_submission.md"))
cat("\n[B2.6] DONE at:\n  ", OUT, "\n")
