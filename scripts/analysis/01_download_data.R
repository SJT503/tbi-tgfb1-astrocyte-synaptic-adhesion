source("scripts/00_utils.R")

args <- parse_args()
dry_run <- is_true(args$dry_run %||% TRUE)
download_supp <- is_true(args$download_supp %||% FALSE)

config <- load_config("config/datasets.yaml")
manifest <- validate_manifest(read_manifest())

ensure_dir(result_dir("logs"))
ensure_dir(result_dir("tables"))

inventory <- data.frame(
  accession = character(),
  modality = character(),
  species = character(),
  role = character(),
  raw_dir = character(),
  manifest_rows = integer(),
  dry_run = logical(),
  stringsAsFactors = FALSE
)

download_messages <- character()

log_message("Download script started. dry_run=", dry_run, ", download_supp=", download_supp)

if (!dry_run) {
  require_or_stop("GEOquery")
}

for (ds in config$datasets) {
  accession <- ds$accession
  accession_dir <- dataset_dir(accession)
  ensure_dir(accession_dir)

  dataset_manifest <- manifest[manifest$dataset == accession, , drop = FALSE]
  dataset_note <- c(
    paste0("Accession: ", accession),
    paste0("Role: ", ds$role %||% "unassigned"),
    paste0("Species: ", ds$species %||% "unknown"),
    paste0("Modality: ", ds$modality %||% "unknown"),
    paste0("Main analysis: ", ds$include_main_analysis %||% FALSE),
    paste0("Notes: ", ds$notes %||% "NA"),
    paste0("Manifest rows: ", nrow(dataset_manifest))
  )
  save_placeholder_message(file.path(accession_dir, "DATASET_INFO.txt"), dataset_note)

  inventory <- rbind(
    inventory,
    data.frame(
      accession = accession,
      modality = ds$modality %||% NA_character_,
      species = ds$species %||% NA_character_,
      role = ds$role %||% NA_character_,
      raw_dir = normalizePath(accession_dir, winslash = "/", mustWork = FALSE),
      manifest_rows = nrow(dataset_manifest),
      dry_run = dry_run,
      stringsAsFactors = FALSE
    )
  )

  if (dry_run) {
    download_messages <- c(download_messages, paste0(accession, ": dry run only, no network request sent."))
    next
  }

  tryCatch(
    {
      GEOquery::getGEO(accession, GSEMatrix = FALSE, destdir = accession_dir)
      if (download_supp) {
        GEOquery::getGEOSuppFiles(accession, baseDir = accession_dir, makeDirectory = FALSE)
      }
      download_messages <- c(download_messages, paste0(accession, ": download completed."))
    },
    error = function(e) {
      download_messages <<- c(download_messages, paste0(accession, ": download failed - ", conditionMessage(e)))
    }
  )
}

write_csv_safe(inventory, result_dir("tables", "dataset_inventory.csv"))
save_placeholder_message(
  result_dir("logs", "download_summary.log"),
  c(
    paste0("dry_run=", dry_run),
    paste0("download_supp=", download_supp),
    download_messages
  )
)

log_message("Download script finished. Inventory written to results/tables/dataset_inventory.csv")
