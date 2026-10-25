source("scripts/00_utils.R")

require_or_stop(c("Seurat", "yaml", "monocle3", "ggplot2"))

params <- load_config("config/analysis_params.yaml")
astro_path <- processed_dir("astrocyte", "astrocyte_object.rds")

if (!file.exists(astro_path)) {
  stop("Missing astrocyte object: ", astro_path)
}

ensure_dir(processed_dir("astrocyte"))
ensure_dir(result_dir("tables"))
ensure_dir(result_dir("logs"))

astro_obj <- readRDS(astro_path)
astro_obj <- attach_inferred_metadata(astro_obj)
state_template_path <- result_dir("tables", "astrocyte_state_template.csv")
if (file.exists(state_template_path)) {
  astro_obj <- apply_astro_state_labels(
    astro_obj,
    utils::read.csv(state_template_path, stringsAsFactors = FALSE)
  )
}
astro_obj <- tryCatch(
  SeuratObject::JoinLayers(astro_obj, assay = "RNA"),
  error = function(e) astro_obj
)
Seurat::DefaultAssay(astro_obj) <- if ("RNA" %in% names(astro_obj@assays)) "RNA" else Seurat::DefaultAssay(astro_obj)

counts_mat <- Seurat::GetAssayData(astro_obj, assay = "RNA", layer = "counts")
gene_meta <- data.frame(
  gene_short_name = rownames(counts_mat),
  row.names = rownames(counts_mat),
  stringsAsFactors = FALSE
)
cell_meta <- astro_obj[[]]

cds <- monocle3::new_cell_data_set(
  expression_data = counts_mat,
  cell_metadata = cell_meta,
  gene_metadata = gene_meta
)
cds <- monocle3::preprocess_cds(cds, num_dim = 30)
cds <- monocle3::reduce_dimension(cds, reduction_method = "UMAP")
cds <- monocle3::cluster_cells(cds, reduction_method = "UMAP")
cds <- monocle3::learn_graph(cds, use_partition = TRUE)

root_cells <- colnames(astro_obj)[astro_obj$astro_state == params$trajectory$root_state_hint]
if (length(root_cells) == 0) {
  root_cells <- colnames(astro_obj)[Seurat::Idents(astro_obj) == levels(Seurat::Idents(astro_obj))[1]]
}
root_cells <- head(root_cells, 500)
cds <- monocle3::order_cells(cds, root_cells = root_cells)

pseudotime_values <- monocle3::pseudotime(cds)
astro_obj$pseudotime <- pseudotime_values[colnames(astro_obj)]
saveRDS(astro_obj, processed_dir("astrocyte", "astrocyte_object_with_pseudotime.rds"))
saveRDS(cds, processed_dir("astrocyte", "astrocyte_monocle3_cds.rds"))

pt_table <- data.frame(
  cell = colnames(astro_obj),
  astro_cluster = as.character(Seurat::Idents(astro_obj)),
  astro_state = astro_obj$astro_state,
  dataset = astro_obj$dataset,
  sample_id = astro_obj$sample_id,
  condition = astro_obj$condition,
  brain_region = astro_obj$brain_region,
  pseudotime = astro_obj$pseudotime,
  stringsAsFactors = FALSE
)

state_summary <- aggregate(
  pseudotime ~ astro_state + condition,
  data = pt_table,
  FUN = function(x) c(n = length(x), mean = mean(x, na.rm = TRUE), median = stats::median(x, na.rm = TRUE))
)
state_summary_expanded <- data.frame(
  astro_state = state_summary$astro_state,
  condition = state_summary$condition,
  n_cells = state_summary$pseudotime[, "n"],
  mean_pseudotime = state_summary$pseudotime[, "mean"],
  median_pseudotime = state_summary$pseudotime[, "median"],
  stringsAsFactors = FALSE
)

expr_mat <- Seurat::GetAssayData(astro_obj, assay = "RNA", layer = "data")
variable_genes <- head(Seurat::VariableFeatures(astro_obj), 1500)
variable_genes <- intersect(variable_genes, rownames(expr_mat))
gene_cor <- lapply(
  variable_genes,
  function(gene) {
    values <- as.numeric(expr_mat[gene, ])
    keep <- is.finite(values) & is.finite(astro_obj$pseudotime)
    if (sum(keep) < 50) {
      return(NULL)
    }
    cor_val <- suppressWarnings(stats::cor(values[keep], astro_obj$pseudotime[keep], method = "spearman"))
    data.frame(gene = gene, spearman_rho = cor_val, stringsAsFactors = FALSE)
  }
)
gene_cor <- do.call(rbind, gene_cor)
gene_cor <- gene_cor[order(-abs(gene_cor$spearman_rho)), , drop = FALSE]

plot_pseudotime <- monocle3::plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE
)
plot_state <- monocle3::plot_cells(
  cds,
  color_cells_by = "astro_state",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE
)

ggplot2::ggsave(
  filename = result_dir("figures", "astrocyte_pseudotime_umap.pdf"),
  plot = plot_pseudotime,
  width = 8,
  height = 6
)
ggplot2::ggsave(
  filename = result_dir("figures", "astrocyte_state_umap_monocle3.pdf"),
  plot = plot_state,
  width = 8,
  height = 6
)

trajectory_status <- data.frame(
  method = params$trajectory$method,
  astro_cells = ncol(astro_obj),
  astro_features = nrow(astro_obj),
  root_state_hint = params$trajectory$root_state_hint,
  status = "completed",
  note = "Monocle3 pseudotime completed on reclustered astrocytes.",
  stringsAsFactors = FALSE
)

write_csv_safe(trajectory_status, result_dir("tables", "trajectory_status.csv"))
write_csv_safe(pt_table, result_dir("tables", "astrocyte_pseudotime_per_cell.csv"))
write_csv_safe(state_summary_expanded, result_dir("tables", "astrocyte_pseudotime_state_summary.csv"))
write_csv_safe(head(gene_cor, 200), result_dir("tables", "astrocyte_pseudotime_top_genes.csv"))
save_placeholder_message(
  result_dir("logs", "trajectory.log"),
  c(
    paste0("Requested method: ", params$trajectory$method),
    paste0("Root state hint: ", params$trajectory$root_state_hint),
    "Monocle3 pseudotime and trajectory plots were generated."
  )
)

log_message("Trajectory analysis completed.")
