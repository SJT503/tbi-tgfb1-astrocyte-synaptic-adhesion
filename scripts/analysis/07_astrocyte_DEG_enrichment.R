source("scripts/00_utils.R")

require_or_stop(c("Seurat", "yaml"))

params <- load_config("config/analysis_params.yaml")
astro_path <- processed_dir("astrocyte", "astrocyte_object.rds")
state_template_path <- result_dir("tables", "astrocyte_state_template.csv")

if (!file.exists(astro_path)) {
  stop("Missing astrocyte object: ", astro_path)
}

ensure_dir(processed_dir("astrocyte"))
ensure_dir(result_dir("tables"))
ensure_dir(result_dir("logs"))

astro_obj <- readRDS(astro_path)
astro_obj <- attach_inferred_metadata(astro_obj)
astro_obj <- tryCatch(
  SeuratObject::JoinLayers(astro_obj, assay = "RNA"),
  error = function(e) astro_obj
)
Seurat::DefaultAssay(astro_obj) <- if ("RNA" %in% names(astro_obj@assays)) "RNA" else Seurat::DefaultAssay(astro_obj)

if (file.exists(state_template_path)) {
  state_template <- utils::read.csv(state_template_path, stringsAsFactors = FALSE)
} else {
  state_template <- data.frame(
    astro_cluster = levels(Seurat::Idents(astro_obj)),
    suggested_state = "manual_review_required",
    stringsAsFactors = FALSE
  )
}
astro_obj <- apply_astro_state_labels(astro_obj, state_template)

deg_results <- tryCatch(
  {
    Seurat::FindAllMarkers(
      astro_obj,
      only.pos = FALSE,
      logfc.threshold = params$deg$logfc_threshold,
      min.pct = params$deg$min_pct
    )
  },
  error = function(e) {
    data.frame(error = conditionMessage(e), stringsAsFactors = FALSE)
  }
)

state_markers <- tryCatch(
  {
    Seurat::FindAllMarkers(
      astro_obj,
      only.pos = TRUE,
      group.by = "astro_state",
      logfc.threshold = params$deg$logfc_threshold,
      min.pct = params$deg$min_pct
    )
  },
  error = function(e) {
    data.frame(error = conditionMessage(e), stringsAsFactors = FALSE)
  }
)

top_state_markers <- if (all(c("cluster", "avg_log2FC") %in% colnames(state_markers))) {
  do.call(
    rbind,
    lapply(
      split(state_markers, state_markers$cluster),
      function(df) head(df[order(-df$avg_log2FC), , drop = FALSE], 20)
    )
  )
} else {
  state_markers
}

state_condition_summary <- as.data.frame.matrix(
  table(
    astro_state = astro_obj$astro_state,
    condition = astro_obj$condition,
    useNA = "ifany"
  )
)
state_condition_summary$astro_state <- rownames(state_condition_summary)
state_condition_summary <- state_condition_summary[
  ,
  c("astro_state", setdiff(colnames(state_condition_summary), "astro_state"))
]

signature_table <- data.frame(
  signature_name = params$astrocyte_states$signature_modules,
  status = "pending_gene_set_input",
  note = "Add custom gene sets under config/gene_sets before module scoring.",
  stringsAsFactors = FALSE
)

write_csv_safe(deg_results, result_dir("tables", "astrocyte_findallmarkers.csv"))
write_csv_safe(top_state_markers, result_dir("tables", "astrocyte_state_top_markers.csv"))
write_csv_safe(state_condition_summary, result_dir("tables", "astrocyte_state_condition_summary.csv"))
write_csv_safe(signature_table, result_dir("tables", "astrocyte_signature_status.csv"))
save_placeholder_message(
  result_dir("logs", "astrocyte_deg_enrichment.log"),
  c(
    paste0("Astrocyte clusters: ", paste(levels(Seurat::Idents(astro_obj)), collapse = ", ")),
    paste0("State template rows: ", nrow(state_template)),
    "Pathway enrichment and module scoring should be run after gene sets are added."
  )
)

log_message("Astrocyte DEG and state marker tables exported.")
