`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) {
    y
  } else {
    x
  }
}

project_root <- normalizePath(".", winslash = "/", mustWork = FALSE)

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

configure_library_paths <- function() {
  local_lib <- file.path(project_root, "local_env", "R", "library")
  ensure_dir(local_lib)

  user_lib <- Sys.getenv("R_LIBS_USER", unset = "")
  lib_paths <- c(local_lib, user_lib, .libPaths())
  lib_paths <- unique(lib_paths[nzchar(lib_paths)])
  .libPaths(lib_paths)

  invisible(.libPaths())
}

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  parsed <- list()
  if (length(args) == 0) {
    return(parsed)
  }

  for (arg in args) {
    if (!startsWith(arg, "--")) {
      next
    }

    body <- sub("^--", "", arg)
    if (grepl("=", body, fixed = TRUE)) {
      parts <- strsplit(body, "=", fixed = TRUE)[[1]]
      key <- parts[1]
      value <- paste(parts[-1], collapse = "=")
      parsed[[key]] <- value
    } else {
      parsed[[body]] <- TRUE
    }
  }

  parsed
}

is_true <- function(x) {
  if (is.logical(x)) {
    return(isTRUE(x))
  }
  tolower(as.character(x) %||% "false") %in% c("true", "t", "1", "yes", "y")
}

load_config <- function(path) {
  yaml::read_yaml(path)
}

save_placeholder_message <- function(path, message) {
  ensure_dir(dirname(path))
  writeLines(message, con = path, useBytes = TRUE)
}

log_message <- function(...) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message("[", timestamp, "] ", paste(..., collapse = ""))
}

write_csv_safe <- function(df, path, row.names = FALSE) {
  ensure_dir(dirname(path))
  utils::write.csv(df, file = path, row.names = row.names, quote = TRUE)
}

read_rds_if_valid <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }

  tryCatch(
    readRDS(path),
    error = function(e) {
      message("Invalid or incomplete RDS detected, ignoring: ", path)
      NULL
    }
  )
}

require_or_stop <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing required packages: ",
      paste(missing, collapse = ", "),
      ". Install them before running this step."
    )
  }
}

read_manifest <- function(path = "config/sample_manifest.csv") {
  if (!file.exists(path)) {
    stop("Missing sample manifest: ", path)
  }
  utils::read.csv(path, stringsAsFactors = FALSE, strip.white = TRUE)
}

validate_manifest <- function(manifest) {
  required_cols <- c(
    "dataset", "sample_id", "species", "modality", "technology", "brain_region",
    "condition", "timepoint", "include_main_analysis", "analysis_tier", "note"
  )
  missing_cols <- setdiff(required_cols, colnames(manifest))
  if (length(missing_cols) > 0) {
    stop("Sample manifest is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  manifest$sample_id <- trimws(manifest$sample_id)
  manifest$species <- trimws(manifest$species)
  manifest$include_main_analysis <- toupper(trimws(as.character(manifest$include_main_analysis)))

  manifest
}

dataset_dir <- function(accession) {
  file.path("raw_data", accession)
}

processed_dir <- function(...) {
  file.path("processed_data", ...)
}

result_dir <- function(...) {
  file.path("results", ...)
}

find_supported_expression_inputs <- function(accession_dir) {
  patterns <- c("\\.h5$", "matrix\\.mtx(\\.gz)?$", "barcodes\\.tsv(\\.gz)?$", "features\\.tsv(\\.gz)?$", "genes\\.tsv(\\.gz)?$")
  found <- unlist(lapply(
    patterns,
    function(pattern) {
      list.files(accession_dir, pattern = pattern, recursive = TRUE, full.names = TRUE)
    }
  ))
  unique(found)
}

detect_species_mt_pattern <- function(species) {
  species_lower <- tolower(species %||% "")
  if (species_lower == "human") {
    "^MT-"
  } else {
    "^mt-"
  }
}

build_cohort_label <- function(brain_region, condition, timepoint) {
  paste(brain_region, condition, timepoint, sep = "_")
}

configure_library_paths()

suppressPackageStartupMessages({
  library(yaml)
})

safe_make_names <- function(x) {
  make.names(x, unique = TRUE)
}

read_gse180862_metadata <- function(accession_dir = dataset_dir("GSE180862")) {
  metadata_files <- list.files(
    accession_dir,
    pattern = "\\.metaData\\.tsv\\.gz$",
    full.names = TRUE
  )

  if (length(metadata_files) == 0) {
    return(data.frame())
  }

  metadata_tables <- lapply(
    metadata_files,
    function(path) {
      prefix <- sub("\\.metaData\\.tsv\\.gz$", "", basename(path))
      meta <- utils::read.delim(gzfile(path), stringsAsFactors = FALSE)
      meta$cell_short <- rownames(meta)
      meta$cell_full <- paste0(prefix, "_", meta$cell_short)
      meta$sample_id <- meta$Sample %||% prefix
      meta$timepoint <- meta$Timepoint %||% NA_character_
      meta$condition_raw <- meta$Condition %||% NA_character_
      meta$condition <- ifelse(
        tolower(meta$condition_raw) %in% c("hn", "tbi", "injury", "injured"),
        "tbi",
        ifelse(
          tolower(meta$condition_raw) %in% c("sham", "control", "ctrl", "intact"),
          "control",
          tolower(meta$condition_raw)
        )
      )
      meta$animal <- meta$Animal %||% NA_character_
      meta$brain_region <- ifelse(
        grepl("hip", meta$Tissue, ignore.case = TRUE),
        "hippocampus",
        ifelse(
          grepl("cortex", meta$Tissue, ignore.case = TRUE),
          "cortex",
          tolower(meta$Tissue %||% NA_character_)
        )
      )
      meta$dataset <- "GSE180862"
      meta
    }
  )

  meta <- do.call(rbind, metadata_tables)
  rownames(meta) <- meta$cell_full
  meta
}

infer_metadata_from_cells <- function(cell_names, dataset_values = NULL) {
  if (is.null(dataset_values)) {
    dataset_values <- ifelse(
      grepl("^GSE180862_", cell_names),
      "GSE180862",
      ifelse(
        grepl("^(Sham|TBI)[0-9]+_", cell_names, ignore.case = TRUE),
        "GSE101901",
        NA_character_
      )
    )
  }

  meta <- data.frame(
    cell = cell_names,
    dataset = dataset_values,
    sample_id = NA_character_,
    condition = NA_character_,
    condition_raw = NA_character_,
    timepoint = NA_character_,
    brain_region = NA_character_,
    animal = NA_character_,
    stringsAsFactors = FALSE
  )
  rownames(meta) <- cell_names

  gse101901_idx <- which(meta$dataset == "GSE101901")
  if (length(gse101901_idx) > 0) {
    prefixes <- sub("_.*$", "", meta$cell[gse101901_idx])
    meta$sample_id[gse101901_idx] <- prefixes
    meta$condition_raw[gse101901_idx] <- ifelse(
      grepl("^sham", prefixes, ignore.case = TRUE),
      "Sham",
      "TBI"
    )
    meta$condition[gse101901_idx] <- ifelse(
      grepl("^sham", prefixes, ignore.case = TRUE),
      "control",
      "tbi"
    )
    meta$timepoint[gse101901_idx] <- "24h"
    meta$brain_region[gse101901_idx] <- "hippocampus"
    meta$animal[gse101901_idx] <- prefixes
  }

  gse180862_idx <- which(meta$dataset == "GSE180862")
  if (length(gse180862_idx) > 0) {
    geo_meta <- read_gse180862_metadata()
    matched <- match(meta$cell[gse180862_idx], rownames(geo_meta))
    found <- !is.na(matched)

    if (any(found)) {
      rows <- geo_meta[matched[found], , drop = FALSE]
      target_idx <- gse180862_idx[found]
      meta$sample_id[target_idx] <- rows$sample_id
      meta$condition[target_idx] <- rows$condition
      meta$condition_raw[target_idx] <- rows$condition_raw
      meta$timepoint[target_idx] <- rows$timepoint
      meta$brain_region[target_idx] <- rows$brain_region
      meta$animal[target_idx] <- rows$animal
    }
  }

  meta
}

attach_inferred_metadata <- function(obj) {
  inferred <- infer_metadata_from_cells(
    cell_names = colnames(obj),
    dataset_values = if ("dataset" %in% colnames(obj[[]])) obj$dataset else NULL
  )

  for (col in colnames(inferred)) {
    if (col == "cell") {
      next
    }
    obj[[col]] <- inferred[[col]]
  }

  obj
}

assign_cluster_celltypes <- function(obj, marker_list, manual_table = NULL) {
  obj <- tryCatch(
    SeuratObject::JoinLayers(obj, assay = "RNA"),
    error = function(e) obj
  )
  Seurat::DefaultAssay(obj) <- if ("RNA" %in% names(obj@assays)) "RNA" else Seurat::DefaultAssay(obj)

  marker_genes <- unique(unlist(marker_list, use.names = FALSE))
  marker_genes <- intersect(marker_genes, rownames(obj))
  avg_expr <- Seurat::AverageExpression(
    obj,
    assays = "RNA",
    features = marker_genes,
    return.seurat = FALSE,
    verbose = FALSE
  )$RNA

  cluster_ids <- colnames(avg_expr)
  cluster_ids <- sub("^g(?=[0-9]+$)", "", cluster_ids, perl = TRUE)
  score_list <- lapply(
    marker_list,
    function(genes) {
      genes <- intersect(genes, rownames(avg_expr))
      if (length(genes) == 0) {
        rep(NA_real_, ncol(avg_expr))
      } else if (length(genes) == 1) {
        as.numeric(avg_expr[genes, , drop = TRUE])
      } else {
        row_subset <- as.matrix(avg_expr[genes, , drop = FALSE])
        colMeans(row_subset)
      }
    }
  )

  score_matrix <- do.call(cbind, score_list)
  colnames(score_matrix) <- names(marker_list)
  rownames(score_matrix) <- cluster_ids

  auto_type <- apply(
    score_matrix,
    1,
    function(x) {
      if (all(is.na(x))) {
        return("manual_review_required")
      }
      names(which.max(x))[1]
    }
  )

  auto_score <- apply(
    score_matrix,
    1,
    function(x) {
      if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
    }
  )

  annotation <- data.frame(
    cluster = cluster_ids,
    suggested_cell_type = auto_type[cluster_ids],
    auto_score = auto_score[cluster_ids],
    rationale = "Auto-assigned from canonical marker average expression.",
    stringsAsFactors = FALSE
  )

  if (!is.null(manual_table) && nrow(manual_table) > 0) {
    manual_table$cluster <- as.character(manual_table$cluster)
    override_idx <- match(annotation$cluster, manual_table$cluster)
    has_override <- !is.na(override_idx) &
      nzchar(manual_table$suggested_cell_type[override_idx]) &
      manual_table$suggested_cell_type[override_idx] != "manual_review_required"
    annotation$suggested_cell_type[has_override] <- manual_table$suggested_cell_type[override_idx[has_override]]
    annotation$rationale[has_override] <- manual_table$rationale[override_idx[has_override]]
  }

  cluster_map <- setNames(annotation$suggested_cell_type, annotation$cluster)
  obj$celltype_auto <- unname(cluster_map[as.character(Seurat::Idents(obj))])

  list(obj = obj, annotation = annotation, score_matrix = score_matrix)
}

apply_astro_state_labels <- function(obj, state_template) {
  state_template$astro_cluster <- as.character(state_template$astro_cluster)
  state_map <- setNames(state_template$suggested_state, state_template$astro_cluster)
  obj$astro_state <- unname(state_map[as.character(Seurat::Idents(obj))])
  obj$astro_state[is.na(obj$astro_state) | !nzchar(obj$astro_state)] <- "manual_review_required"
  obj
}
