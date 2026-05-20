#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(SCENT)
})

parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected positional argument: ", key)
    }
    name <- sub("^--", "", key)
    if (name %in% c("force")) {
      out[[name]] <- TRUE
      i <- i + 1
    } else {
      if (i == length(args)) {
        stop("Missing value for ", key)
      }
      out[[name]] <- args[[i + 1]]
      i <- i + 2
    }
  }
  out
}

required_arg <- function(args, name) {
  value <- args[[name]]
  if (is.null(value) || !nzchar(value)) {
    stop("Missing required argument --", name)
  }
  value
}

read_gz_lines <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con), add = TRUE)
  readLines(con, warn = FALSE)
}

read_gz_table <- function(path, ...) {
  con <- gzfile(path, "rt")
  on.exit(close(con), add = TRUE)
  read.delim(con, check.names = FALSE, stringsAsFactors = FALSE, ...)
}

read_sparse_matrix <- function(path, row_names, col_names) {
  con <- gzfile(path, "rt")
  on.exit(close(con), add = TRUE)
  mat <- Matrix::readMM(con)
  mat <- as(mat, "CsparseMatrix")
  if (nrow(mat) != length(row_names)) {
    stop("Matrix row count mismatch for ", path, ": got ", nrow(mat), " rows but ", length(row_names), " row names.")
  }
  if (ncol(mat) != length(col_names)) {
    stop("Matrix column count mismatch for ", path, ": got ", ncol(mat), " columns but ", length(col_names), " cell names.")
  }
  rownames(mat) <- row_names
  colnames(mat) <- col_names
  mat
}

slugify <- function(value) {
  slug <- gsub("[^A-Za-z0-9]+", "_", value)
  slug <- gsub("^_+|_+$", "", slug)
  if (!nzchar(slug)) "celltype" else slug
}

current_rss_mb <- function() {
  status_path <- "/proc/self/status"
  if (!file.exists(status_path)) {
    return(NA_real_)
  }
  status <- readLines(status_path, warn = FALSE)
  vmrss <- grep("^VmRSS:", status, value = TRUE)
  if (length(vmrss) == 0) {
    return(NA_real_)
  }
  as.numeric(gsub("[^0-9]", "", vmrss[[1]])) / 1024
}

log_step <- function(...) {
  message(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), " | ", ...)
}

script_start <- Sys.time()
args <- parse_args(commandArgs(trailingOnly = TRUE))

prepared_dir <- required_arg(args, "prepared-dir")
results_dir <- required_arg(args, "results-dir")
logs_dir <- required_arg(args, "logs-dir")
batch <- as.integer(required_arg(args, "batch"))
celltype <- required_arg(args, "celltype")
cores <- as.integer(ifelse(is.null(args[["cores"]]), "6", args[["cores"]]))
regr <- ifelse(is.null(args[["regr"]]), "poisson", args[["regr"]])
bin <- as.logical(ifelse(is.null(args[["bin"]]), "TRUE", args[["bin"]]))
force <- isTRUE(args[["force"]])
downsample_cells <- as.integer(ifelse(is.null(args[["downsample-cells"]]), "0", args[["downsample-cells"]]))
downsample_seed <- as.integer(ifelse(is.null(args[["downsample-seed"]]), "1", args[["downsample-seed"]]))

batch_label <- sprintf("batch_%03d", batch)
batch_dir <- file.path(prepared_dir, "batches", batch_label)
peak_info_path <- file.path(prepared_dir, "peak_info_batches", paste0("peak_info_", batch_label, ".tsv.gz"))
celltype_slug <- slugify(celltype)
celltype_outdir <- file.path(results_dir, celltype_slug)
dir.create(celltype_outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(logs_dir, showWarnings = FALSE, recursive = TRUE)

outfile <- file.path(celltype_outdir, paste0(batch_label, ".csv"))
failfile <- file.path(celltype_outdir, paste0(batch_label, ".failed.txt"))

if (file.exists(outfile) && file.info(outfile)$size > 0 && !force) {
  message("Output already exists for ", celltype, " ", batch_label, "; skipping: ", outfile)
  quit(save = "no", status = 0)
}
if (file.exists(failfile) && force) {
  unlink(failfile)
}

load_start <- Sys.time()
log_step("Loading prepared SCENT inputs for ", celltype, " ", batch_label)
cells <- read_gz_lines(file.path(prepared_dir, "cells.tsv.gz"))
genes <- read_gz_lines(file.path(prepared_dir, "genes.tsv.gz"))
peaks <- read_gz_lines(file.path(batch_dir, "atac_peaks.tsv.gz"))

metadata <- read_gz_table(file.path(prepared_dir, "metadata.tsv.gz"))
if (!identical(metadata$cell_id, cells)) {
  stop("metadata.tsv.gz cell_id order does not match cells.tsv.gz.")
}
metadata$cell <- metadata$cell_id
rownames(metadata) <- metadata$cell_id

peak_info <- read_gz_table(peak_info_path)
peak_info <- peak_info[, c("Gene", "Peak")]

rna_matrix <- read_sparse_matrix(file.path(prepared_dir, "rna_selected_genes.mtx.gz"), genes, cells)
atac_matrix <- read_sparse_matrix(file.path(batch_dir, "atac.mtx.gz"), peaks, cells)
bin_for_scent <- bin
if (isTRUE(bin)) {
  log_step("Pre-binarizing sparse ATAC matrix before SCENT to avoid SCENT zero-peak binning failures.")
  atac_matrix@x[atac_matrix@x > 0] <- 1
  bin_for_scent <- FALSE
}
log_step("Loaded matrices in ", round(as.numeric(difftime(Sys.time(), load_start, units = "mins")), 2),
         " minutes; current RSS MB=", round(current_rss_mb(), 1))

missing_genes <- setdiff(unique(peak_info$Gene), rownames(rna_matrix))
missing_peaks <- setdiff(unique(peak_info$Peak), rownames(atac_matrix))
if (length(missing_genes) > 0) {
  stop("peak_info contains genes missing from RNA matrix; first examples: ", paste(head(missing_genes, 10), collapse = ", "))
}
if (length(missing_peaks) > 0) {
  stop("peak_info contains peaks missing from ATAC matrix; first examples: ", paste(head(missing_peaks, 10), collapse = ", "))
}
if (!celltype %in% metadata$celltype) {
  stop("Requested cell type is not present in metadata: ", celltype)
}

covariates <- c("n_counts", "log1p_n_fragment", "tsse", "frac_mito", "frac_dup")
missing_covariates <- setdiff(covariates, colnames(metadata))
if (length(missing_covariates) > 0) {
  stop("Metadata is missing covariates: ", paste(missing_covariates, collapse = ", "))
}

celltype_idx <- which(metadata$celltype == celltype)
if (length(celltype_idx) == 0) {
  stop("No cells found for requested cell type: ", celltype)
}

if (downsample_cells > 0 && length(celltype_idx) > downsample_cells) {
  set.seed(downsample_seed)
  keep_idx <- sort(sample(celltype_idx, size = downsample_cells, replace = FALSE))
  log_step("Downsampling cell type ", celltype, " from ", length(celltype_idx),
           " to ", length(keep_idx), " cells with seed=", downsample_seed)
  metadata <- metadata[keep_idx, , drop = FALSE]
  selected_cells <- metadata$cell_id
  rna_matrix <- rna_matrix[, selected_cells, drop = FALSE]
  atac_matrix <- atac_matrix[, selected_cells, drop = FALSE]
} else {
  log_step("Using all ", length(celltype_idx), " cells for cell type ", celltype)
}

log_step("Creating SCENT object with ", nrow(peak_info), " gene-peak pairs, ",
         nrow(rna_matrix), " RNA genes, ", nrow(atac_matrix), " ATAC bins, ",
         ncol(rna_matrix), " selected cells.")

status <- tryCatch({
  object_start <- Sys.time()
  scent_obj <- CreateSCENTObj(
    rna = rna_matrix,
    atac = atac_matrix,
    meta.data = metadata,
    peak.info = peak_info,
    covariates = covariates,
    celltypes = "celltype"
  )
  log_step("Created SCENT object in ", round(as.numeric(difftime(Sys.time(), object_start, units = "mins")), 2),
           " minutes; current RSS MB=", round(current_rss_mb(), 1))

  options(warn = -1)
  start_time <- Sys.time()
  log_step("Starting SCENT_algorithm; cores=", cores, ", regr=", regr,
           ", requested_bin=", bin, ", passed_bin=", bin_for_scent)
  scent_obj <- SCENT_algorithm(
    object = scent_obj,
    celltype = celltype,
    ncores = cores,
    regr = regr,
    bin = bin_for_scent
  )
  elapsed <- round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 2)

  result <- scent_obj@SCENT.result
  write.csv(result, outfile, row.names = FALSE)
  log_step("Finished SCENT_algorithm for ", celltype, " ", batch_label, " in ", elapsed,
           " minutes; result rows=", nrow(result), "; current RSS MB=", round(current_rss_mb(), 1))
  log_step("Total R script runtime minutes=", round(as.numeric(difftime(Sys.time(), script_start, units = "mins")), 2),
           "; output=", outfile)
  0
}, error = function(e) {
  msg <- paste0("SCENT failed for celltype=", celltype, " batch=", batch, "\n", conditionMessage(e), "\n")
  writeLines(msg, failfile)
  log_step(msg)
  log_step("Total R script runtime before failure minutes=", round(as.numeric(difftime(Sys.time(), script_start, units = "mins")), 2),
           "; current RSS MB=", round(current_rss_mb(), 1))
  1
})

quit(save = "no", status = status)
