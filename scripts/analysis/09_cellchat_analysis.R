source("scripts/00_utils.R")

require_or_stop(c("Seurat", "yaml", "CellChat"))

params <- load_config("config/analysis_params.yaml")
integrated_path <- processed_dir("integrated", "mouse_integrated_annotated.rds")

if (!file.exists(integrated_path)) {
  stop("Missing annotated integrated object: ", integrated_path)
}

ensure_dir(processed_dir("communication"))
ensure_dir(result_dir("tables"))
ensure_dir(result_dir("logs"))

obj <- readRDS(integrated_path)
obj <- attach_inferred_metadata(obj)

annotation_path <- result_dir("tables", "cluster_annotation_template.csv")
annotation_table <- if (file.exists(annotation_path)) {
  utils::read.csv(annotation_path, stringsAsFactors = FALSE)
} else {
  NULL
}

annotation_result <- assign_cluster_celltypes(
  obj = obj,
  marker_list = params$annotation$major_celltype_markers,
  manual_table = annotation_table
)
obj <- annotation_result$obj

focus_types <- unique(c(
  params$communication$focus_sender,
  params$communication$focus_receivers
))
obj$condition_binary <- ifelse(
  obj$condition %in% c("control", "sham", "intact"),
  "control",
  ifelse(obj$condition %in% c("tbi", "injury", "injured"), "tbi", obj$condition)
)

keep <- obj$celltype_auto %in% focus_types & obj$condition_binary %in% c("control", "tbi")
obj <- subset(obj, cells = colnames(obj)[keep])
obj <- tryCatch(
  SeuratObject::JoinLayers(obj, assay = "RNA"),
  error = function(e) obj
)
Seurat::DefaultAssay(obj) <- if ("RNA" %in% names(obj@assays)) "RNA" else Seurat::DefaultAssay(obj)

avg_expr <- Seurat::AggregateExpression(
  obj,
  assays = "RNA",
  group.by = c("condition_binary", "celltype_auto"),
  return.seurat = FALSE,
  verbose = FALSE
)$RNA

get_complex_score <- function(gene_string, expr_vector) {
  genes <- strsplit(gene_string, "_", fixed = TRUE)[[1]]
  genes <- genes[genes %in% names(expr_vector)]
  if (length(genes) == 0) {
    return(NA_real_)
  }
  mean(expr_vector[genes], na.rm = TRUE)
}

score_interactions_for_condition <- function(condition_label) {
  sender_label <- paste(condition_label, params$communication$focus_sender, sep = "_")
  if (!sender_label %in% colnames(avg_expr)) {
    return(data.frame())
  }

  sender_expr <- avg_expr[, sender_label, drop = TRUE]
  interaction_db <- CellChat::CellChatDB.mouse$interaction
  results <- list()
  idx <- 0L

  for (receiver in params$communication$focus_receivers) {
    receiver_label <- paste(condition_label, receiver, sep = "_")
    if (!receiver_label %in% colnames(avg_expr)) {
      next
    }
    receiver_expr <- avg_expr[, receiver_label, drop = TRUE]

    for (i in seq_len(nrow(interaction_db))) {
      ligand_score <- get_complex_score(interaction_db$ligand[i], sender_expr)
      receptor_score <- get_complex_score(interaction_db$receptor[i], receiver_expr)
      if (!is.finite(ligand_score) || !is.finite(receptor_score)) {
        next
      }

      idx <- idx + 1L
      results[[idx]] <- data.frame(
        condition = condition_label,
        source = params$communication$focus_sender,
        target = receiver,
        interaction_name = interaction_db$interaction_name[i],
        pathway_name = interaction_db$pathway_name[i],
        ligand = interaction_db$ligand[i],
        receptor = interaction_db$receptor[i],
        ligand_expr = ligand_score,
        receptor_expr = receptor_score,
        score = ligand_score * receptor_score,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(results) == 0) {
    return(data.frame())
  }

  out <- do.call(rbind, results)
  out$priority_match <- out$pathway_name %in% params$communication$priority_pathways
  out[order(-out$score), , drop = FALSE]
}

control_df <- score_interactions_for_condition("control")
tbi_df <- score_interactions_for_condition("tbi")
all_interactions <- rbind(control_df, tbi_df)

comparison_df <- merge(
  control_df[, c("source", "target", "interaction_name", "pathway_name", "score")],
  tbi_df[, c("source", "target", "interaction_name", "pathway_name", "score")],
  by = c("source", "target", "interaction_name", "pathway_name"),
  all = TRUE,
  suffixes = c("_control", "_tbi")
)
comparison_df$score_control[is.na(comparison_df$score_control)] <- 0
comparison_df$score_tbi[is.na(comparison_df$score_tbi)] <- 0
comparison_df$delta_tbi_minus_control <- comparison_df$score_tbi - comparison_df$score_control
comparison_df$priority_match <- comparison_df$pathway_name %in% params$communication$priority_pathways
comparison_df <- comparison_df[order(-comparison_df$delta_tbi_minus_control), , drop = FALSE]

write_csv_safe(
  head(control_df, 200),
  result_dir("tables", "cellchat_like_control_top200.csv")
)
write_csv_safe(
  head(tbi_df, 200),
  result_dir("tables", "cellchat_like_tbi_top200.csv")
)
write_csv_safe(
  head(comparison_df, 300),
  result_dir("tables", "cellchat_astro_sender_priority.csv")
)
write_csv_safe(all_interactions, result_dir("tables", "cellchat_all_interactions.csv"))

saveRDS(comparison_df, processed_dir("communication", "cellchat_like_comparison.rds"))

communication_plan <- data.frame(
  method = paste0(params$communication$method, "_DB_scored"),
  focus_sender = params$communication$focus_sender,
  focus_receivers = paste(params$communication$focus_receivers, collapse = ";"),
  priority_pathways = paste(params$communication$priority_pathways, collapse = ";"),
  status = "completed",
  note = "Used CellChat mouse ligand-receptor database with average-expression scoring for robust first-pass communication ranking.",
  stringsAsFactors = FALSE
)

write_csv_safe(communication_plan, result_dir("tables", "communication_plan.csv"))
save_placeholder_message(
  result_dir("logs", "cellchat.log"),
  c(
    paste0("Planned method: ", params$communication$method),
    paste0("Focus sender: ", params$communication$focus_sender),
    paste0("Retained cells: ", ncol(obj)),
    "Communication ranking was generated with CellChatDB-backed ligand-receptor scoring."
  )
)

log_message("CellChat-style communication ranking completed.")
