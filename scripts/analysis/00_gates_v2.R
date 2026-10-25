## 00_gates_v2.R — TBI 严重度梯度项目 v2 门禁函数库
## PROJECT_CONSTITUTION_v2.md 第四节定义
## 所有 scripts/200+ 必须 source 此文件，禁止直接 readRDS()

PROJ_ROOT <- "D:/2026-04-24/new-chat-2/TBI_astrocyte_scst_project"
V2_DATA   <- file.path(PROJ_ROOT, "processed_data", "v2")

# ── 唯一 RDS 注册表 ──────────────────────────────────────────────────
V2_RDS_REGISTRY <- list(
  severity_integrated = file.path(V2_DATA, "severity_integrated.rds"),
  human_ast           = file.path(V2_DATA, "human_gse209552_ast.rds"),
  spatial             = file.path(V2_DATA, "spatial_gse226208.rds")
)

# ── 读取门控函数 ─────────────────────────────────────────────────────
load_v2_rds <- function(name, expect_min_cells = 500, check_purity = TRUE,
                        purity_ratio = 3.0) {
  path <- V2_RDS_REGISTRY[[name]]
  if (is.null(path)) stop("GATE FAIL: '", name, "' 不在 V2_RDS_REGISTRY")
  if (!file.exists(path)) stop("GATE FAIL: 文件不存在 → ", path)

  message("Loading: ", path)
  seu <- readRDS(path)

  # 检验 ① 细胞数
  if (ncol(seu) < expect_min_cells)
    stop("GATE FAIL: 细胞数 ", ncol(seu), " < ", expect_min_cells)

  # 检验 ② 必要元数据列
  required_cols <- c("severity", "dataset", "condition")
  missing_cols  <- setdiff(required_cols, colnames(seu@meta.data))
  if (length(missing_cols) > 0)
    stop("GATE FAIL: 缺少元数据列: ", paste(missing_cols, collapse = ", "))

  # 检验 ③ 纯度（仅星形胶质 RDS）
  if (check_purity && name %in% c("severity_integrated")) {
    expr <- tryCatch(
      Seurat::GetAssayData(seu, layer = "data"),
      error = function(e) Seurat::GetAssayData(seu, slot = "data")
    )
    genes_present <- rownames(expr)
    if ("Aqp4" %in% genes_present && "Snap25" %in% genes_present) {
      ratio <- mean(expr["Aqp4", ]) / (mean(expr["Snap25", ]) + 1e-6)
      if (ratio < purity_ratio)
        stop("PURITY FAIL: Aqp4/Snap25 = ", round(ratio, 2),
             " < ", purity_ratio, " — 不是纯星形胶质 RDS")
      message("Purity OK: Aqp4/Snap25 = ", round(ratio, 2))
    }
  }

  message("Gate PASS: ", name, " (", ncol(seu), " cells, ",
          nlevels(factor(seu@meta.data$severity)), " severity levels)")
  return(seu)
}

# ── 严重度标准标签（全局统一）───────────────────────────────────────
SEVERITY_LEVELS  <- c("sham", "mild", "moderate", "severe")
SEVERITY_COLORS  <- c(sham     = "#4DAECD",
                       mild     = "#8BC4D9",
                       moderate = "#E07B54",
                       severe   = "#B03A2E")

DATASET_SEVERITY <- c(
  GSE247339  = "mild",
  GSE101901  = "mild",
  GSE290150  = "moderate",
  GSE180862  = "moderate",
  GSE226207  = "moderate",
  GSE271769  = "severe",
  GSE277487  = "severe"
)

# ── 亚态标准标签（现代多亚态框架）──────────────────────────────────
SUBSTATE_LEVELS <- c("homeostatic", "inflammatory", "IFN_responsive",
                     "reactive_gliotic", "DAA", "border_associated",
                     "metabolic_remodeled")
SUBSTATE_MARKERS <- list(
  homeostatic       = c("Aldoc", "Glul", "Aqp4", "Slc1a3"),
  inflammatory      = c("Lcn2", "Steap4", "Cxcl10", "Serpina3n"),
  IFN_responsive    = c("Ifit1", "Ifit2", "Ifit3", "Mx1", "Oasl2"),
  reactive_gliotic  = c("Gfap", "Vim", "Nestin", "Cd44"),
  DAA               = c("Gfap", "Apoe", "Ctsb", "Ctsd"),
  border_associated = c("Igfbp2", "Gpc3", "Ptn"),
  metabolic_remodeled = c("Hk2", "Pfkfb3", "Ldha")
)

message("gates_v2.R loaded — V2 project constitution active")
