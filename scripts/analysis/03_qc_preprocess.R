source("scripts/00_utils.R")

args <- parse_args()
target_dataset <- args$dataset %||% NA_character_

params <- load_config("config/analysis_params.yaml")
manifest <- validate_manifest(read_manifest())
manifest$include_main_analysis <- toupper(manifest$include_main_analysis) == "TRUE"

ensure_dir(processed_dir("qc"))
ensure_dir(result_dir("qc_reports"))
ensure_dir(result_dir("tables"))
ensure_dir(result_dir("logs"))

read_gse101901_txt <- function(file_path) {
  require_or_stop(c("Matrix"))

  con <- gzfile(file_path, open = "rt")
  on.exit(close(con), add = TRUE)

  barcode_line <- readLines(con, n = 1)
  barcodes <- strsplit(barcode_line, "\t", fixed = TRUE)[[1]]

  i_idx <- integer()
  j_idx <- integer()
  x_vals <- numeric()
  gene_ids <- character()
  row_counter <- 0L

  repeat {
    lines <- readLines(con, n = 1000)
    if (length(lines) == 0) {
      break
    }

    for (line in lines) {
      fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
      if (length(fields) != length(barcodes) + 1) {
        stop("GSE101901 row length does not match barcode header structure.")
      }

      row_counter <- row_counter + 1L
      gene_ids[row_counter] <- safe_make_names(fields[1])
      values <- as.numeric(fields[-1])
      nz <- which(values != 0)

      if (length(nz) > 0) {
        i_idx <- c(i_idx, rep.int(row_counter, length(nz)))
        j_idx <- c(j_idx, nz)
        x_vals <- c(x_vals, values[nz])
      }
    }
  }

  Matrix::sparseMatrix(
    i = i_idx,
    j = j_idx,
    x = x_vals,
    dims = c(length(gene_ids), length(barcodes)),
    dimnames = list(gene_ids, barcodes)
  )
}

read_geo_mtx_bundle <- function(dataset_id, accession_dir) {
  require_or_stop(c("Matrix"))

  matrix_files <- list.files(accession_dir, pattern = "digital_expression\\.mtx\\.gz$", full.names = TRUE)
  if (length(matrix_files) == 0) {
    stop("No digital_expression.mtx.gz files found for ", dataset_id)
  }

  bundle_objects <- list()
  for (matrix_file in matrix_files) {
    prefix <- sub("\\.digital_expression\\.mtx\\.gz$", "", basename(matrix_file))
    barcode_file <- file.path(accession_dir, paste0(prefix, ".barcodes.tsv.gz"))
    feature_file <- file.path(accession_dir, paste0(prefix, ".features.tsv.gz"))

    if (!file.exists(barcode_file) || !file.exists(feature_file)) {
      stop("Missing matching barcode/features files for prefix: ", prefix)
    }

    counts <- Matrix::readMM(gzfile(matrix_file))
    barcodes <- readLines(gzfile(barcode_file))
    features <- read.delim(gzfile(feature_file), header = FALSE, stringsAsFactors = FALSE)

    if (ncol(counts) == nrow(features) && nrow(counts) != nrow(features)) {
      counts <- Matrix::t(counts)
    }

    if (nrow(features) < nrow(counts)) {
      stop("Feature rows are fewer than matrix rows for prefix: ", prefix)
    }

    rownames(counts) <- safe_make_names(features[[1]][seq_len(nrow(counts))])
    colnames(counts) <- paste0(prefix, "_", safe_make_names(barcodes))
    bundle_objects[[prefix]] <- counts
  }

  all_genes <- unique(unlist(lapply(bundle_objects, rownames)))

  bundle_objects <- lapply(
    bundle_objects,
    function(mat) {
      if (length(setdiff(all_genes, rownames(mat))) > 0) {
        missing_genes <- setdiff(all_genes, rownames(mat))
        zero_block <- Matrix::Matrix(
          0,
          nrow = length(missing_genes),
          ncol = ncol(mat),
          sparse = TRUE,
          dimnames = list(missing_genes, colnames(mat))
        )
        mat <- rbind(mat, zero_block)
      }
      mat[all_genes, , drop = FALSE]
    }
  )

  do.call(cbind, bundle_objects)
}

load_counts_for_dataset <- function(dataset_id, accession_dir) {
  txt_files <- list.files(accession_dir, pattern = "digital_expression\\.txt\\.gz$", full.names = TRUE)
  if (length(txt_files) > 0) {
    return(read_gse101901_txt(txt_files[1]))
  }

  matrix_files <- list.files(accession_dir, pattern = "digital_expression\\.mtx\\.gz$", full.names = TRUE)
  if (length(matrix_files) > 0) {
    return(read_geo_mtx_bundle(dataset_id, accession_dir))
  }

  stop("No supported count matrix files found for dataset: ", dataset_id)
}

run_qc_for_dataset <- function(dataset_id, dataset_manifest, params) {
  accession_dir <- dataset_dir(dataset_id)
  supported_inputs <- find_supported_expression_inputs(accession_dir)

  status <- data.frame(
    dataset = dataset_id,
    status = "pending",
    n_manifest_rows = nrow(dataset_manifest),
    n_input_files = length(supported_inputs),
    output_rds = NA_character_,
    note = NA_character_,
    stringsAsFactors = FALSE
  )

  if (length(supported_inputs) == 0 && length(list.files(accession_dir, pattern = "digital_expression\\.txt\\.gz$", full.names = TRUE)) == 0) {
    status$status <- "skipped"
    status$note <- "No supported expression files detected under raw_data."
    return(status)
  }

  if (!requireNamespace("Seurat", quietly = TRUE)) {
    status$status <- "skipped"
    status$note <- "Seurat is not installed."
    return(status)
  }

  tryCatch(
    {
      counts <- load_counts_for_dataset(dataset_id, accession_dir)

      obj <- Seurat::CreateSeuratObject(
        counts = counts,
        project = dataset_id,
        min.cells = 3,
        min.features = 0
      )

      species <- unique(dataset_manifest$species)
      species <- species[!is.na(species) & nzchar(species)]
      mt_pattern <- detect_species_mt_pattern(species[1] %||% "mouse")

      obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = mt_pattern)

      max_mt <- if (tolower(species[1] %||% "mouse") == "human") {
        params$qc$max_percent_mt_human
      } else {
        params$qc$max_percent_mt_mouse
      }

      obj <- subset(
        obj,
        subset =
          nFeature_RNA >= params$qc$min_features &
          nFeature_RNA <= params$qc$max_features &
          percent.mt <= max_mt
      )

      output_rds <- processed_dir("qc", paste0(dataset_id, "_qc.rds"))
      saveRDS(obj, output_rds)

      qc_summary <- data.frame(
        dataset = dataset_id,
        cells_after_qc = ncol(obj),
        genes_after_qc = nrow(obj),
        max_percent_mt = max_mt,
        stringsAsFactors = FALSE
      )
      write_csv_safe(qc_summary, result_dir("tables", paste0(dataset_id, "_qc_summary.csv")))

      status$status <- "completed"
      status$output_rds <- normalizePath(output_rds, winslash = "/", mustWork = FALSE)
      status$note <- "QC completed using GEO count matrix loader."
      status
    },
    error = function(e) {
      status$status <- "failed"
      status$note <- conditionMessage(e)
      status
    }
  )
}

dataset_ids <- unique(manifest$dataset[manifest$include_main_analysis])
if (!is.na(target_dataset)) {
  dataset_ids <- intersect(dataset_ids, target_dataset)
}

qc_status <- if (length(dataset_ids) == 0) {
  data.frame(
    dataset = character(),
    status = character(),
    n_manifest_rows = integer(),
    n_input_files = integer(),
    output_rds = character(),
    note = character(),
    stringsAsFactors = FALSE
  )
} else {
  do.call(
    rbind,
    lapply(
      dataset_ids,
      function(dataset_id) {
        dataset_manifest <- manifest[manifest$dataset == dataset_id, , drop = FALSE]
        run_qc_for_dataset(dataset_id, dataset_manifest, params)
      }
    )
  )
}

write_csv_safe(qc_status, result_dir("tables", "qc_status_summary.csv"))
save_placeholder_message(
  result_dir("logs", "qc_preprocess.log"),
  if (nrow(qc_status) == 0) {
    "No datasets matched the QC run criteria."
  } else {
    apply(qc_status, 1, function(row) paste(row, collapse = " | "))
  }
)

log_message("QC preprocessing finished for ", nrow(qc_status), " dataset(s).")
