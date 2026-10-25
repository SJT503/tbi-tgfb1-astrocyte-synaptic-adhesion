source("scripts/00_utils.R")

require_or_stop(c("Seurat", "yaml"))

args <- parse_args()
params <- load_config("config/analysis_params.yaml")
integrated_path <- processed_dir("integrated", "mouse_integrated_clustered.rds")

if (!file.exists(integrated_path)) {
  stop("Missing integrated object: ", integrated_path)
}

ensure_dir(processed_dir("integrated"))
ensure_dir(result_dir("tables"))
ensure_dir(result_dir("logs"))

obj <- readRDS(integrated_path)
obj <- attach_inferred_metadata(obj)
Seurat::DefaultAssay(obj) <- if ("RNA" %in% names(obj@assays)) "RNA" else Seurat::DefaultAssay(obj)
obj <- tryCatch(
  SeuratObject::JoinLayers(obj, assay = "RNA"),
  error = function(e) obj
)

markers <- params$annotation$major_celltype_markers
celltype_summary <- data.frame(
  cell_type = names(markers),
  marker_genes = vapply(markers, paste, collapse = ";", character(1)),
  stringsAsFactors = FALSE
)

cluster_ids <- levels(Seurat::Idents(obj))
manual_template_path <- result_dir("tables", "cluster_annotation_template.csv")
manual_template <- if (file.exists(manual_template_path)) {
  utils::read.csv(manual_template_path, stringsAsFactors = FALSE)
} else {
  NULL
}

annotation_result <- assign_cluster_celltypes(obj, markers, manual_template)
obj <- annotation_result$obj
annotation_table <- annotation_result$annotation

marker_path <- result_dir("tables", "cluster_markers_top10.csv")
marker_matrix_path <- result_dir("tables", "canonical_marker_cluster_average.csv")
cluster_score_path <- result_dir("tables", "cluster_celltype_scores.csv")
force_markers <- is_true(args$force_markers)
marker_results <- if (file.exists(marker_path) && !force_markers) {
  "reused_existing"
} else {
  tryCatch(
    {
      markers <- Seurat::FindAllMarkers(
        object = obj,
        only.pos = TRUE,
        min.pct = 0.25,
        logfc.threshold = 0.25
      )
      markers <- markers[order(markers$cluster, -markers$avg_log2FC), , drop = FALSE]
      top_markers <- do.call(
        rbind,
        lapply(
          split(markers, markers$cluster),
          function(df) head(df, 10)
        )
      )
      write_csv_safe(top_markers, marker_path)
      "completed"
    },
    error = function(e) {
      save_placeholder_message(
        result_dir("logs", "annotation_marker_error.log"),
        conditionMessage(e)
      )
      "failed"
    }
  )
}

marker_genes <- unique(unlist(markers, use.names = FALSE))
marker_genes <- intersect(marker_genes, rownames(obj))
if (length(marker_genes) > 0) {
  avg_expr <- Seurat::AverageExpression(
    obj,
    assays = "RNA",
    features = marker_genes,
    return.seurat = FALSE,
    verbose = FALSE
  )
  avg_expr_df <- as.data.frame(avg_expr$RNA)
  avg_expr_df$gene <- rownames(avg_expr_df)
  avg_expr_df <- avg_expr_df[, c("gene", setdiff(colnames(avg_expr_df), "gene"))]
  write_csv_safe(avg_expr_df, marker_matrix_path)
}

output_path <- processed_dir("integrated", "mouse_integrated_annotated.rds")
saveRDS(obj, output_path)
write_csv_safe(celltype_summary, result_dir("tables", "canonical_marker_reference.csv"))
write_csv_safe(annotation_table, result_dir("tables", "cluster_annotation_template.csv"))
cluster_score_df <- as.data.frame(annotation_result$score_matrix)
cluster_score_df$cluster <- rownames(cluster_score_df)
cluster_score_df <- cluster_score_df[, c("cluster", setdiff(colnames(cluster_score_df), "cluster"))]
write_csv_safe(cluster_score_df, cluster_score_path)
save_placeholder_message(
  result_dir("logs", "annotation.log"),
  c(
    "Annotation template generated.",
    paste0("Marker export status: ", marker_results),
    "Cluster-level canonical marker scores were exported for fast review."
  )
)

log_message("Cell annotation exported with inferred metadata and cluster-level auto labels.")
