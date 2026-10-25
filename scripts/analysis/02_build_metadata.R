source("scripts/00_utils.R")

config <- load_config("config/datasets.yaml")
manifest <- validate_manifest(read_manifest())

ensure_dir(result_dir("tables"))
ensure_dir(result_dir("logs"))

manifest$cohort_label <- build_cohort_label(
  manifest$brain_region,
  manifest$condition,
  manifest$timepoint
)
manifest$is_placeholder <- grepl("placeholder", manifest$sample_id, ignore.case = TRUE)
manifest$include_main_analysis <- toupper(manifest$include_main_analysis) == "TRUE"

dataset_roles <- vapply(
  config$datasets,
  function(ds) ds$role %||% NA_character_,
  character(1)
)
names(dataset_roles) <- vapply(config$datasets, function(ds) ds$accession, character(1))
manifest$dataset_role <- dataset_roles[manifest$dataset]

main_manifest <- manifest[manifest$include_main_analysis, , drop = FALSE]

cohort_summary <- aggregate(
  sample_id ~ dataset + dataset_role + species + brain_region + condition + timepoint + analysis_tier,
  data = manifest,
  FUN = length
)
colnames(cohort_summary)[colnames(cohort_summary) == "sample_id"] <- "n_rows"

placeholder_rows <- manifest[manifest$is_placeholder, , drop = FALSE]
placeholder_message <- if (nrow(placeholder_rows) > 0) {
  c(
    "Placeholder sample IDs detected.",
    "Replace these with exact GEO sample or file identifiers before running downstream analyses.",
    paste(placeholder_rows$dataset, placeholder_rows$sample_id, sep = ": ")
  )
} else {
  "No placeholder sample IDs detected."
}

write_csv_safe(manifest, result_dir("tables", "sample_manifest_checked.csv"))
write_csv_safe(main_manifest, result_dir("tables", "sample_manifest_main_analysis.csv"))
write_csv_safe(cohort_summary, "config/cohort_summary.csv")
save_placeholder_message(result_dir("logs", "metadata_summary.log"), placeholder_message)

log_message("Metadata build finished. Placeholder rows: ", nrow(placeholder_rows))
log_message("Main-analysis rows exported to results/tables/sample_manifest_main_analysis.csv")
