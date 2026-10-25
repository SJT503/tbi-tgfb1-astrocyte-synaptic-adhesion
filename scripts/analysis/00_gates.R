# ==============================================================================
# 00_gates.R — 全项目通用门禁函数库
# ==============================================================================
# 用法：每个分析脚本顶部加 source("00_gates.R")，然后调用：
#   load_astrocyte_rds()  → 自动门禁验证返回 Seurat 对象
#   load_neuron_rds("excitatory") / load_neuron_rds("inhibitory")
#   load_microglia_rds()
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
})

# ── 定稿RDS对照表（焊死） ────────────────────────────────────────────────────
DEFINITIVE_RDS <- list(
  astrocyte  = list(
    file     = "astrocyte_v2_unified.rds",
    n_cells  = 774,
    n_clusters = 3,
    purity_marker = "Aqp4",
    contam_marker = "Snap25",
    min_ratio = 3,
    required_groups = c("Sham","CCI_Scr","CCI_TGFbetaR1_KO"),
    expected_labels = c("Homeostatic-A","Homeostatic-B(Apoe+)","Reactive(Id3+)")
  ),
  excitatory = list(
    file     = "excitatory_v2_unified.rds",
    n_cells  = 18556,
    n_clusters = 15,
    purity_marker = "Snap25",
    contam_marker = "Aqp4",
    min_ratio = 3,
    required_groups = c("Sham","CCI_Scr","CCI_TGFbetaR1_KO")
  ),
  inhibitory = list(
    file     = "inhibitory_v2_unified.rds",
    n_cells  = 3677,
    n_clusters = 9,
    purity_marker = "Snap25",
    contam_marker = "Aqp4",
    min_ratio = 3,
    required_groups = c("Sham","CCI_Scr","CCI_TGFbetaR1_KO")
  ),
  microglia  = list(
    file     = "microglia_v2_unified.rds",
    n_cells  = 2126,
    n_clusters = 6,
    purity_marker = "Hexb",
    contam_marker = "Snap25",
    min_ratio = 3,
    required_groups = c("Sham","CCI_Scr","CCI_TGFbetaR1_KO")
  )
)

SEURAT_DIR <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project/wet_lab_data/results/seurat"

# ── 核心函数：带门禁的RDS加载器 ──────────────────────────────────────────────
load_definitive_rds <- function(celltype) {
  if (!celltype %in% names(DEFINITIVE_RDS))
    stop(sprintf("[GATE FAIL] 未知细胞类型: %s。支持的类型: %s",
                 celltype, paste(names(DEFINITIVE_RDS), collapse=", ")))

  spec <- DEFINITIVE_RDS[[celltype]]
  rds_path <- file.path(SEURAT_DIR, spec$file)

  cat(sprintf("[LOAD] 读取 %s（%s）...\n", celltype, spec$file))
  if (!file.exists(rds_path))
    stop(sprintf("[GATE FAIL] RDS 缺失: %s\n本项目唯一定稿 RDS 由 PROJECT_CONSTITUTION.md 规定，请勿手动指定其他文件。", rds_path))

  seu <- readRDS(rds_path)

  # ── 门禁 1：细胞数 ──
  if (!is.null(spec$n_cells) && ncol(seu) != spec$n_cells)
    stop(sprintf("[GATE FAIL #1] 细胞数错误！期望 %d，实际 %d。RDS可能被替换。",
                 spec$n_cells, ncol(seu)))

  # ── 门禁 2：分群数 ──
  n_k <- length(unique(seu$seurat_clusters))
  if (n_k != spec$n_clusters)
    stop(sprintf("[GATE FAIL #2] 分群数错误！期望 %d 群，实际 %d 群。",
                 spec$n_clusters, n_k))

  # ── 门禁 3：纯度检验 ──
  expr <- tryCatch(GetAssayData(seu, layer="data", assay="RNA"),
                   error = function(e) GetAssayData(seu, slot="data", assay="RNA"))
  if (spec$purity_marker %in% rownames(expr) && spec$contam_marker %in% rownames(expr)) {
    ratio <- mean(expr[spec$purity_marker,]) / max(mean(expr[spec$contam_marker,]), 1e-9)
    if (ratio < spec$min_ratio)
      stop(sprintf("[GATE FAIL #3] 纯度失败！%s/%s = %.2f（需≥%d）",
                   spec$purity_marker, spec$contam_marker, ratio, spec$min_ratio))
    cat(sprintf("  纯度: %s/%s = %.2f ✓\n",
                spec$purity_marker, spec$contam_marker, ratio))
  }

  # ── 门禁 4：分组完整性 ──
  group_col <- intersect(c("stim","group"), colnames(seu@meta.data))[1]
  if (is.na(group_col))
    stop("[GATE FAIL #4] metadata 缺 stim/group 列")
  groups <- table(seu@meta.data[[group_col]])
  missing <- setdiff(spec$required_groups, names(groups))
  if (length(missing) > 0)
    stop(sprintf("[GATE FAIL #4] 缺分组: %s", paste(missing, collapse=", ")))

  # ── 门禁 5：标注完整性（若期望） ──
  if (!is.null(spec$expected_labels)) {
    if (!"cluster_label" %in% colnames(seu@meta.data))
      stop("[GATE FAIL #5] 缺 cluster_label！必须用 63_annotate_clusters.R 标注后的 RDS")
    actual_labels <- unique(seu$cluster_label)
    missing_labels <- setdiff(spec$expected_labels, actual_labels)
    if (length(missing_labels) > 0)
      stop(sprintf("[GATE FAIL #5] 缺标注: %s",
                   paste(missing_labels, collapse=", ")))
    cat(sprintf("  标注: %s ✓\n", paste(spec$expected_labels, collapse=", ")))
  }

  cat(sprintf("[GATE PASS] %d cells, %d clusters, %d groups\n",
              ncol(seu), n_k, length(groups)))

  invisible(seu)
}

# 便捷函数
load_astrocyte_rds  <- function() load_definitive_rds("astrocyte")
load_neuron_rds     <- function(subtype = c("excitatory","inhibitory")) {
  load_definitive_rds(match.arg(subtype))
}
load_microglia_rds  <- function() load_definitive_rds("microglia")

# ── DEG 合理性检验 ───────────────────────────────────────────────────────────
sanity_check_deg <- function(deg_df, gene, expected_direction,
                              gene_col = "gene", fc_col = "avg_log2FC") {
  row <- deg_df[deg_df[[gene_col]] == gene, ]
  if (nrow(row) == 0) {
    cat(sprintf("  [SANITY] %s 不在 DEG 列表（NA）\n", gene))
    return(invisible(NA))
  }
  actual <- sign(row[[fc_col]][1])
  pass <- actual == expected_direction
  cat(sprintf("  [SANITY] %-12s log2FC=%+.3f %s (期望方向=%+d)\n",
              gene, row[[fc_col]][1], ifelse(pass,"✓","⚠"), expected_direction))
  invisible(pass)
}

# ── 出口审计函数 ─────────────────────────────────────────────────────────────
audit_outputs <- function(figout_dir, script_name = "Unknown") {
  cat("\n========== 出口审计 ==========\n")
  cat(sprintf("脚本: %s\n", script_name))
  cat(sprintf("输出目录: %s\n", figout_dir))
  if (dir.exists(figout_dir)) {
    n_pdf <- length(list.files(figout_dir, pattern="\\.pdf$"))
    n_png <- length(list.files(figout_dir, pattern="\\.png$"))
    n_csv <- length(list.files(figout_dir, pattern="\\.csv$"))
    cat(sprintf("产出: %d PDF, %d PNG, %d CSV\n", n_pdf, n_png, n_csv))
  } else {
    cat("  输出目录不存在！\n")
  }
  cat("===============================\n")
}

cat("[00_gates.R] 门禁函数库已加载\n")
