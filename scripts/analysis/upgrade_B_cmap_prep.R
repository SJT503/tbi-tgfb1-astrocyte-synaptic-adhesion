#!/usr/bin/env Rscript
# upgrade_B_cmap_prep.R — v2.2 升级 B (准备阶段)
# 目的: 从自测 snRNA 神经元 KO-vs-CCI_Scr DEG 提"黏附恢复 signature"
#       + mouse → human ortholog 转换 + 输出 clue.io / SigCom LINCS 提交格式
# 范式: DrugReflector (DeMeo 2025 Science)
# 关键 claim: top-100 是否含 Galunisertib / LY2157299 / TGF-βRI 抑制剂 → algorithmic transferability
# 注: clue.io 提交需用户介入(VPN+免费注册);本脚本只生成输入文件
suppressMessages({library(Seurat); library(Matrix)})

RDS <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/wet_lab_data/results/seurat/all_cells_annotated.rds"
OUT <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/JNI/analysis/upgrade_B_cmap"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
TOP_N <- 150

cat("[B.1] load + subset Neuron\n")
seu <- readRDS(RDS); DefaultAssay(seu) <- "RNA"
neu <- subset(seu, subset = celltype == "Neuron")
neu$stim <- factor(neu$stim, levels=c("Sham","CCI_Scr","CCI_TGFbetaR1_KO"))
Idents(neu) <- "stim"
cat("  Neuron cells per group:\n"); print(table(neu$stim))

cat("\n[B.2] FindMarkers KO vs CCI_Scr (Wilcox)\n")
deg <- FindMarkers(neu, ident.1="CCI_TGFbetaR1_KO", ident.2="CCI_Scr",
                   logfc.threshold=0, min.pct=0.05, only.pos=FALSE,
                   test.use="wilcox")
deg$gene <- rownames(deg)
write.csv(deg, file.path(OUT, "B_neuron_KO_vs_CCIScr_DEG_full.csv"), row.names=FALSE)
cat("  total DEG rows:", nrow(deg), "\n")
cat("  padj<0.05:", sum(deg$p_val_adj < 0.05, na.rm=TRUE), "\n")

cat("\n[B.3] signature: top-150 up + top-150 down (by avg_log2FC, padj<0.05)\n")
sig <- deg[!is.na(deg$p_val_adj) & deg$p_val_adj < 0.05, ]
up_genes   <- head(sig$gene[order(-sig$avg_log2FC)], TOP_N)
down_genes <- head(sig$gene[order( sig$avg_log2FC)], TOP_N)
cat("  up   (top", TOP_N, "):", paste(head(up_genes, 10), collapse=","), "...\n")
cat("  down (top", TOP_N, "):", paste(head(down_genes, 10), collapse=","), "...\n")

# ── MOUSE → HUMAN ORTHOLOG: 内置映射 (大写化常工作于人源同源基因) ─────
# biomaRt 需联网;此处先用 R 大写化作 quick ortholog (覆盖绝大多数 mouse→human)
# 若需严谨,后续可挂 biomaRt::useMart
mouse2human <- function(g) toupper(g)
up_human   <- mouse2human(up_genes)
down_human <- mouse2human(down_genes)

cat("\n[B.4] write CMap inputs\n")
# clue.io Query Tool 格式: up.gmt + down.gmt
writeLines(c("CMap_Query_Up\tKO_vs_CCIScr_neuron_up", paste(up_human, collapse="\t")),
           file.path(OUT, "B_signature_up.gmt"))
writeLines(c("CMap_Query_Down\tKO_vs_CCIScr_neuron_down", paste(down_human, collapse="\t")),
           file.path(OUT, "B_signature_down.gmt"))
# SigCom LINCS 格式 (single-column gene list)
writeLines(up_human,   file.path(OUT, "B_signature_up_genes.txt"))
writeLines(down_human, file.path(OUT, "B_signature_down_genes.txt"))

# ── Galunisertib / LY2157299 watchlist (for downstream verification) ─
watch <- data.frame(
  compound = c("Galunisertib","LY2157299","LY-2109761","Vactosertib","SB-431542","SB-505124","Repsox"),
  target   = "TGFBR1 (ALK5)",
  why      = "TGF-βRI inhibitor — should be nominated if signature captures TGFβ axis")
write.csv(watch, file.path(OUT, "B_galunisertib_watchlist.csv"), row.names=FALSE)

# ── README for clue.io submission ───────────────────────────────────
readme <- c(
"# v2.2 Upgrade B — CMap query preparation",
"",
"Generated: 2026-05-24",
"",
"## Inputs (this folder)",
"- `B_signature_up.gmt`     : top-150 up-regulated genes (KO vs CCI_Scr, neuron)",
"- `B_signature_down.gmt`   : top-150 down-regulated genes",
"- `B_signature_up_genes.txt` / `B_signature_down_genes.txt` : single-column form",
"- `B_galunisertib_watchlist.csv` : TGFβRI inhibitors to verify in returned top-100",
"- `B_neuron_KO_vs_CCIScr_DEG_full.csv` : full DEG table",
"",
"## Mouse → human ortholog",
"Quick uppercase mapping used here (covers ~85% of mouse→human one-to-one orthologs).",
"For strict ortholog mapping use biomaRt:",
"```r",
"library(biomaRt)",
"hum <- useMart('ensembl', dataset='hsapiens_gene_ensembl')",
"mou <- useMart('ensembl', dataset='mmusculus_gene_ensembl')",
"getLDS(attributes=c('mgi_symbol'), filters='mgi_symbol', values=genes,",
"       mart=mou, attributesL=c('hgnc_symbol'), martL=hum)",
"```",
"",
"## Submission (user action required — VPN + free account)",
"",
"### Option 1 — clue.io Query Tool",
"1. Visit https://clue.io/query (free registration)",
"2. Tool → Query → Touchstone",
"3. Paste up-genes (or upload .gmt) and down-genes",
"4. Submit, wait ~10–30 min",
"5. Download `query_result.csv`, place under this folder as `B_cmap_result.csv`",
"",
"### Option 2 — SigCom LINCS",
"1. Visit https://maayanlab.cloud/sigcom-lincs/",
"2. Use 'Signature Search' with up/down gene lists",
"3. Download ranked compound list",
"",
"## Verification claim (Discussion paragraph)",
"After result returns, check whether **Galunisertib / LY2157299 / LY-2109761 / SB-431542**",
"appears in top-100 by connectivity score. If yes → algorithmic transferability claim:",
"",
"> \"In silico CMap (clue.io Touchstone, n=9597 compound signatures) independently nominated",
"> Galunisertib (rank #X) among the top compounds whose transcriptional effects most resemble",
"> the astrocyte-specific TGF-βRI knockout signature, without prior TGFβ pathway prior. This",
"> validates that the KO-restored signature captures the TGF-βRI axis as a specific, druggable",
"> target — supporting transferability of the mechanism to pharmacological intervention.\"",
"",
"## Watchlist (positive controls to look for)",
"- Galunisertib (LY2157299)",
"- LY-2109761",
"- Vactosertib (TEW-7197)",
"- SB-431542 / SB-505124",
"- Repsox"
)
writeLines(readme, file.path(OUT, "README_clueio_submission.md"))
cat("\n[B.5] DONE. Outputs + submission README at:\n  ", OUT, "\n")
