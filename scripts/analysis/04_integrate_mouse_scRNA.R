source("scripts/00_utils.R")

require_or_stop(c("Seurat", "yaml"))

args <- parse_args()
target_datasets <- args$datasets %||% ""
target_datasets <- strsplit(target_datasets, ",", fixed = TRUE)[[1]]
target_datasets <- trimws(target_datasets)
target_datasets <- target_datasets[nzchar(target_datasets)]

params <- load_config("config/analysis_params.yaml")
manifest <- validate_manifest(read_manifest())
manifest$include_main_analysis <- toupper(manifest$include_main_analysis) == "TRUE"

candidate_datasets <- unique(
  manifest$dataset[
    manifest$include_main_analysis &
      manifest$species == "mouse" &
      manifest$modality == "scRNA"
  ]
)

if (length(target_datasets) > 0) {
  candidate_datasets <- intersect(candidate_datasets, target_datasets)
}

if (length(candidate_datasets) == 0) {
  stop("No eligible mouse scRNA datasets found for integration.")
}

ensure_dir(processed_dir("integrated"))
ensure_dir(result_dir("tables"))
ensure_dir(result_dir("logs"))

preprocessed_path <- processed_dir("integrated", "mouse_objects_preprocessed.rds")
features_path <- processed_dir("integrated", "integration_features.rds")
anchors_path <- processed_dir("integrated", "integration_anchors.rds")
integrated_path <- processed_dir("integrated", "mouse_integrated.rds")

objects <- list()
integration_input <- data.frame(
  dataset = character(),
  file = character(),
  cells = integer(),
  features = integer(),
  stringsAsFactors = FALSE
)

for (dataset_id in candidate_datasets) {
  qc_path <- processed_dir("qc", paste0(dataset_id, "_qc.rds"))
  if (!file.exists(qc_path)) {
    log_message("Skipping ", dataset_id, " because QC object is missing: ", qc_path)
    next
  }

  obj <- readRDS(qc_path)
  obj$dataset <- dataset_id
  if (!"sample_id" %in% colnames(obj[[]])) {
    obj$sample_id <- paste0(dataset_id, "_sample")
  }
  objects[[dataset_id]] <- obj

  integration_input <- rbind(
    integration_input,
    data.frame(
      dataset = dataset_id,
      file = normalizePath(qc_path, winslash = "/", mustWork = FALSE),
      cells = ncol(obj),
      features = nrow(obj),
      stringsAsFactors = FALSE
    )
  )
}

if (length(objects) < 2) {
  stop("At least two QC-filtered mouse scRNA objects are required for integration.")
}

preprocessed_cached <- read_rds_if_valid(preprocessed_path)
if (!is.null(preprocessed_cached)) {
  objects <- preprocessed_cached
} else {
  objects <- lapply(
    objects,
    function(obj) {
      Seurat::DefaultAssay(obj) <- "RNA"
      obj <- Seurat::NormalizeData(obj, verbose = FALSE)
      obj <- Seurat::FindVariableFeatures(obj, nfeatures = params$integration$nfeatures, verbose = FALSE)
      obj
    }
  )
  saveRDS(objects, preprocessed_path)
}

features_cached <- read_rds_if_valid(features_path)
if (!is.null(features_cached)) {
  features_to_integrate <- features_cached
} else {
  variable_feature_sets <- lapply(objects, Seurat::VariableFeatures)
  features_to_integrate <- Reduce(intersect, variable_feature_sets)

  if (length(features_to_integrate) < 500) {
    features_to_integrate <- Reduce(intersect, lapply(objects, rownames))
    features_to_integrate <- head(features_to_integrate, params$integration$nfeatures)
  }

  if (length(features_to_integrate) == 0) {
    stop("No shared integration features were identified across the selected objects.")
  }
  saveRDS(features_to_integrate, features_path)
}

anchors_cached <- read_rds_if_valid(anchors_path)
if (!is.null(anchors_cached)) {
  anchors <- anchors_cached
} else {
  anchors <- Seurat::FindIntegrationAnchors(
    object.list = objects,
    anchor.features = features_to_integrate,
    dims = 1:30
  )
  saveRDS(anchors, anchors_path)
}

integrated_cached <- read_rds_if_valid(integrated_path)
if (!is.null(integrated_cached)) {
  integrated <- integrated_cached
} else {
  integrated <- Seurat::IntegrateData(
    anchorset = anchors,
    dims = 1:30
  )
  saveRDS(integrated, integrated_path)
}

scaled_path <- processed_dir("integrated", "mouse_integrated_scaled.rds")
pca_path <- processed_dir("integrated", "mouse_integrated_pca.rds")
umap_path <- processed_dir("integrated", "mouse_integrated_umap.rds")
clustered_path <- processed_dir("integrated", "mouse_integrated_clustered.rds")

scaled_cached <- read_rds_if_valid(scaled_path)
if (!is.null(scaled_cached)) {
  integrated <- scaled_cached
} else {
  Seurat::DefaultAssay(integrated) <- "integrated"
  integrated <- Seurat::ScaleData(integrated, verbose = FALSE)
  saveRDS(integrated, scaled_path)
}

pca_cached <- read_rds_if_valid(pca_path)
if (!is.null(pca_cached)) {
  integrated <- pca_cached
} else {
  integrated <- Seurat::RunPCA(integrated, verbose = FALSE)
  saveRDS(integrated, pca_path)
}

umap_cached <- read_rds_if_valid(umap_path)
if (!is.null(umap_cached)) {
  integrated <- umap_cached
} else {
  integrated <- Seurat::RunUMAP(integrated, dims = 1:30, verbose = FALSE)
  saveRDS(integrated, umap_path)
}

clustered_cached <- read_rds_if_valid(clustered_path)
if (!is.null(clustered_cached)) {
  integrated <- clustered_cached
} else {
  integrated <- Seurat::FindNeighbors(integrated, dims = 1:30, verbose = FALSE)
  integrated <- Seurat::FindClusters(integrated, resolution = params$integration$resolution_all_cells, verbose = FALSE)
  saveRDS(integrated, clustered_path)
}

write_csv_safe(integration_input, result_dir("tables", "integration_input_summary.csv"))
save_placeholder_message(
  result_dir("logs", "integration.log"),
  c(
    paste0("Integrated datasets: ", paste(names(objects), collapse = ", ")),
    paste0("Preprocessed objects: ", normalizePath(preprocessed_path, winslash = "/", mustWork = FALSE)),
    paste0("Integration features: ", normalizePath(features_path, winslash = "/", mustWork = FALSE)),
    paste0("Anchors: ", normalizePath(anchors_path, winslash = "/", mustWork = FALSE)),
    paste0("Integrated raw: ", normalizePath(integrated_path, winslash = "/", mustWork = FALSE)),
    paste0("Scaled: ", normalizePath(scaled_path, winslash = "/", mustWork = FALSE)),
    paste0("PCA: ", normalizePath(pca_path, winslash = "/", mustWork = FALSE)),
    paste0("UMAP: ", normalizePath(umap_path, winslash = "/", mustWork = FALSE)),
    paste0("Clustered output: ", normalizePath(clustered_path, winslash = "/", mustWork = FALSE))
  )
)

log_message("Mouse integration finished: ", clustered_path)
