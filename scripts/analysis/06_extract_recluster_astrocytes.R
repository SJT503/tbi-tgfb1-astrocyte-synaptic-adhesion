source("scripts/00_utils.R")

require_or_stop(c("Seurat", "yaml"))

params <- load_config("config/analysis_params.yaml")
annotated_path <- processed_dir("integrated", "mouse_integrated_annotated.rds")
annotation_template_path <- result_dir("tables", "cluster_annotation_template.csv")

if (!file.exists(annotated_path)) {
  stop("Missing annotated integrated object: ", annotated_path)
}

if (!file.exists(annotation_template_path)) {
  stop("Missing cluster annotation template: ", annotation_template_path)
}

ensure_dir(processed_dir("astrocyte"))
ensure_dir(result_dir("tables"))
ensure_dir(result_dir("logs"))

obj <- readRDS(annotated_path)
obj <- attach_inferred_metadata(obj)
annotation_table <- utils::read.csv(annotation_template_path, stringsAsFactors = FALSE)

astro_clusters <- annotation_table$cluster[tolower(annotation_table$suggested_cell_type) == "astrocyte"]
if (length(astro_clusters) == 0) {
  stop("No clusters are labeled as astrocyte in cluster_annotation_template.csv.")
}

astro_obj <- subset(obj, idents = astro_clusters)

astro_obj <- tryCatch(
  SeuratObject::JoinLayers(astro_obj, assay = "RNA"),
  error = function(e) astro_obj
)
Seurat::DefaultAssay(astro_obj) <- if ("RNA" %in% names(astro_obj@assays)) "RNA" else Seurat::DefaultAssay(astro_obj)
astro_obj <- Seurat::NormalizeData(astro_obj, verbose = FALSE)
astro_obj <- Seurat::FindVariableFeatures(astro_obj, nfeatures = params$integration$nfeatures, verbose = FALSE)
astro_obj <- Seurat::ScaleData(astro_obj, verbose = FALSE)
astro_obj <- Seurat::RunPCA(astro_obj, verbose = FALSE)
astro_obj <- Seurat::RunUMAP(astro_obj, dims = 1:20, verbose = FALSE)
astro_obj <- Seurat::FindNeighbors(astro_obj, dims = 1:20, verbose = FALSE)
astro_obj <- Seurat::FindClusters(astro_obj, resolution = params$integration$resolution_astrocyte, verbose = FALSE)

astro_state_template <- data.frame(
  astro_cluster = levels(Seurat::Idents(astro_obj)),
  suggested_state = "manual_review_required",
  planned_state_pool = paste(params$astrocyte_states$planned_labels, collapse = ";"),
  stringsAsFactors = FALSE
)

existing_state_path <- result_dir("tables", "astrocyte_state_template.csv")
if (file.exists(existing_state_path)) {
  existing_state <- utils::read.csv(existing_state_path, stringsAsFactors = FALSE)
  match_idx <- match(astro_state_template$astro_cluster, existing_state$astro_cluster)
  keep_idx <- !is.na(match_idx)
  astro_state_template$suggested_state[keep_idx] <- existing_state$suggested_state[match_idx[keep_idx]]
}

astro_obj <- apply_astro_state_labels(astro_obj, astro_state_template)

output_path <- processed_dir("astrocyte", "astrocyte_object.rds")
saveRDS(astro_obj, output_path)
write_csv_safe(astro_state_template, result_dir("tables", "astrocyte_state_template.csv"))
save_placeholder_message(
  result_dir("logs", "astrocyte_recluster.log"),
  c(
    paste0("Astrocyte clusters used: ", paste(astro_clusters, collapse = ", ")),
    "Review astrocyte_state_template.csv and assign biological labels after marker analysis."
  )
)

log_message("Astrocyte reclustering finished with state labels: ", output_path)
