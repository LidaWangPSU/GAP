#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# GitHub-shareable TATE / CMap drug-repurposing pipeline
# ------------------------------------------------------------------------------
# This script scores LINCS/CMap perturbations against TWAS gene-level effects.
# It is adapted from an internal pipeline by replacing hard-coded cluster paths
# with command-line arguments and adding basic input validation.
#
# Example:
# Rscript run_tate_pipeline_github_shared.R \
#   --twas-summary /path/to/TWAS_summary_BIOVU_union+100snps-5e-04-lrt_01.txt \
#   --model-name BIOVU_union+100snps-5e-04-lrt_01 \
#   --tissue Whole_Blood \
#   --cell ALL \
#   --pip 0.9 \
#   --cmap-gctx /path/to/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx \
#   --pert-info /path/to/GSE92742_Broad_LINCS_pert_info.txt.gz \
#   --validation-drugs /path/to/cmap.novel.summary.csv \
#   --ensembl-entrez /path/to/ensembl_entrez_id.txt \
#   --processed-cmap-template /path/to/processed_cmap_data/{cell}_20_{chunk}_processed_data.txt \
#   --output-dir /path/to/output
#
# Expected input formats:
#   TWAS summary: columns ID, CHR, BEST.GWAS.ID, TWAS.Z, TWAS.P, tissue
#   ID conversion: columns ensembl_id, entrez_id
#   Validation drugs: column drug, matching pert_iname in pert-info
#   Pert-info: columns pert_id, pert_iname
#   Processed CMap chunks: first column gene/entrez ID, remaining columns perturbations
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(cmapR)
})

option_list <- list(
  make_option("--twas-summary", type = "character", help = "Path to TWAS summary file."),
  make_option("--model-name", type = "character", default = NULL,
              help = "Model name used in output file name. Defaults to TWAS file basename."),
  make_option("--tissue", type = "character", default = NULL,
              help = "Tissue/model name to analyze. If omitted, use --tissue-index."),
  make_option("--tissue-index", type = "integer", default = NULL,
              help = "1-based index into unique TWAS tissue values. Used only if --tissue is omitted."),
  make_option("--cell", type = "character", default = "ALL",
              help = "CMap cell line name used in processed CMap chunk filenames [default: %default]."),
  make_option("--pip", type = "double", default = 0.9,
              help = "Minimum PIP threshold from loci clustering [default: %default]."),
  make_option("--z-threshold", type = "double", default = 0,
              help = "Set CMap values with abs(value) <= threshold to zero [default: %default]."),
  make_option("--n-bootstrap", type = "integer", default = 5000,
              help = "Number of bootstrap/noise replicates [default: %default]."),
  make_option("--chunks", type = "character", default = "1:20",
              help = "CMap chunk indices, e.g. '1:20' or '1,2,3' [default: %default]."),
  make_option("--cmap-gctx", type = "character", help = "Path to LINCS/CMap Level 5 .gctx file."),
  make_option("--pert-info", type = "character", help = "Path to LINCS/CMap perturbation info file."),
  make_option("--validation-drugs", type = "character", help = "CSV with a column named 'drug'."),
  make_option("--ensembl-entrez", type = "character", help = "Gene ID mapping file with ensembl_id and entrez_id columns."),
  make_option("--processed-cmap-template", type = "character",
              help = "Template for processed CMap chunks. Use {cell} and {chunk}, e.g. /dir/{cell}_20_{chunk}_processed_data.txt."),
  make_option("--output-dir", type = "character", default = ".",
              help = "Output directory [default: %default]."),
  make_option("--output-prefix", type = "character", default = "tate_cmap_bootstrap",
              help = "Output filename prefix [default: %default]."),
  make_option("--seed", type = "integer", default = 1,
              help = "Base random seed [default: %default].")
)

opt <- parse_args(OptionParser(option_list = option_list))

stop_if_missing <- function(x, label) {
  if (is.null(x) || is.na(x) || !nzchar(x)) {
    stop(label, " is required.", call. = FALSE)
  }
}

check_file <- function(path, label) {
  stop_if_missing(path, label)
  if (!file.exists(path)) {
    stop(label, " does not exist: ", path, call. = FALSE)
  }
}

parse_chunks <- function(x) {
  x <- gsub("\\s+", "", x)
  if (grepl(":", x)) {
    bounds <- as.integer(strsplit(x, ":")[[1]])
    if (length(bounds) != 2 || any(is.na(bounds))) stop("Invalid --chunks value: ", x, call. = FALSE)
    return(seq(bounds[1], bounds[2]))
  }
  out <- as.integer(strsplit(x, ",")[[1]])
  if (any(is.na(out))) stop("Invalid --chunks value: ", x, call. = FALSE)
  out
}

fill_template <- function(template, cell, chunk) {
  out <- gsub("\\{cell\\}", cell, template)
  out <- gsub("\\{chunk\\}", as.character(chunk), out)
  out
}

# Local replacement for the original sourced GAP.main.function.R::find_loci().
# It groups TWAS signals into loci by chromosome and genomic distance, then keeps
# the top TWAS signal per locus and assigns a simple locus-level PIP proxy if a
# PIP column is not already available.
find_loci_twas <- function(gwas, window = 1e6) {
  required <- c("chr", "pos", "pval")
  missing <- setdiff(required, names(gwas))
  if (length(missing) > 0) stop("find_loci_twas missing columns: ", paste(missing, collapse = ", "), call. = FALSE)

  gwas <- gwas[!is.na(gwas$chr) & !is.na(gwas$pos) & !is.na(gwas$pval), ]
  if (nrow(gwas) == 0) return(gwas)

  gwas <- gwas[order(gwas$chr, gwas$pos, gwas$pval), ]
  gwas$locus <- NA_integer_
  locus_id <- 0L

  for (chr_i in sort(unique(gwas$chr))) {
    idx <- which(gwas$chr == chr_i)
    chr_pos <- gwas$pos[idx]
    current_start <- NA_real_
    current_end <- NA_real_
    for (ii in seq_along(idx)) {
      row_i <- idx[ii]
      pos_i <- chr_pos[ii]
      if (is.na(current_start) || pos_i > current_end) {
        locus_id <- locus_id + 1L
        current_start <- pos_i - window
        current_end <- pos_i + window
      } else {
        current_end <- max(current_end, pos_i + window)
      }
      gwas$locus[row_i] <- locus_id
    }
  }

  gwas$top <- FALSE
  for (locus_i in unique(gwas$locus)) {
    idx <- which(gwas$locus == locus_i)
    top_idx <- idx[which.min(gwas$pval[idx])]
    gwas$top[top_idx] <- TRUE
  }

  if (!"pip" %in% names(gwas)) {
    # Approximate relative posterior weight within each locus from TWAS Z-scores.
    # This is not a fine-mapping PIP; replace with externally computed PIP when available.
    if ("TWAS.Z" %in% names(gwas)) {
      gwas$pip <- NA_real_
      for (locus_i in unique(gwas$locus)) {
        idx <- which(gwas$locus == locus_i)
        z2 <- gwas$TWAS.Z[idx]^2
        weights <- exp(0.5 * (z2 - max(z2, na.rm = TRUE)))
        gwas$pip[idx] <- weights / sum(weights, na.rm = TRUE)
      }
    } else {
      gwas$pip <- as.numeric(gwas$top)
    }
  }
  gwas
}

check_file(opt$twas_summary, "--twas-summary")
check_file(opt$cmap_gctx, "--cmap-gctx")
check_file(opt$pert_info, "--pert-info")
check_file(opt$validation_drugs, "--validation-drugs")
check_file(opt$ensembl_entrez, "--ensembl-entrez")
stop_if_missing(opt$processed_cmap_template, "--processed-cmap-template")

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(opt$output_dir)) stop("Could not create --output-dir: ", opt$output_dir, call. = FALSE)

if (is.null(opt$model_name) || !nzchar(opt$model_name)) {
  opt$model_name <- tools::file_path_sans_ext(basename(opt$twas_summary))
}

message("Reading CMap row metadata...")
row_meta <- read_gctx_meta(opt$cmap_gctx, dim = "row")

message("Reading validation drug list and perturbation metadata...")
validation <- read.csv(opt$validation_drugs, stringsAsFactors = FALSE)
if (!"drug" %in% names(validation)) stop("--validation-drugs must contain a 'drug' column.", call. = FALSE)

drug <- as.data.frame(fread(opt$pert_info, header = TRUE))
required_drug_cols <- c("pert_id", "pert_iname")
missing_drug_cols <- setdiff(required_drug_cols, names(drug))
if (length(missing_drug_cols) > 0) {
  stop("--pert-info missing columns: ", paste(missing_drug_cols, collapse = ", "), call. = FALSE)
}
drug_sub <- drug[drug$pert_iname %in% validation$drug, ]
if (nrow(drug_sub) == 0) warning("No perturbations in --pert-info matched validation$drug.")

message("Reading TWAS summary and gene ID map...")
df <- as.data.frame(fread(opt$twas_summary))
required_twas_cols <- c("ID", "CHR", "BEST.GWAS.ID", "TWAS.Z", "TWAS.P", "tissue")
missing_twas_cols <- setdiff(required_twas_cols, names(df))
if (length(missing_twas_cols) > 0) {
  stop("--twas-summary missing columns: ", paste(missing_twas_cols, collapse = ", "), call. = FALSE)
}

ref <- as.data.frame(fread(opt$ensembl_entrez, header = TRUE))
required_ref_cols <- c("ensembl_id", "entrez_id")
missing_ref_cols <- setdiff(required_ref_cols, names(ref))
if (length(missing_ref_cols) > 0) {
  stop("--ensembl-entrez missing columns: ", paste(missing_ref_cols, collapse = ", "), call. = FALSE)
}
ref <- ref[ref$entrez_id %in% row_meta$id, ]

df$entrez_id <- ref$entrez_id[match(df$ID, ref$ensembl_id)]
df <- df[!is.na(df$entrez_id), ]
df$chr <- df$CHR
df$pos <- suppressWarnings(as.numeric(vapply(strsplit(df$BEST.GWAS.ID, "_"), function(x) x[2], character(1))))
df$pval <- as.numeric(df$TWAS.P)

if (nrow(df) == 0) stop("No TWAS genes remained after Ensembl-to-Entrez mapping.", call. = FALSE)

tissues <- unique(df$tissue)
if (!is.null(opt$tissue) && nzchar(opt$tissue)) {
  tissue_to_run <- opt$tissue
  if (!tissue_to_run %in% tissues) {
    stop("Requested --tissue not found in TWAS summary: ", tissue_to_run,
         "\nAvailable tissues: ", paste(tissues, collapse = ", "), call. = FALSE)
  }
} else {
  if (is.null(opt$tissue_index)) stop("Provide either --tissue or --tissue-index.", call. = FALSE)
  if (opt$tissue_index < 1 || opt$tissue_index > length(tissues)) {
    stop("--tissue-index out of range. There are ", length(tissues), " tissues.", call. = FALSE)
  }
  tissue_to_run <- tissues[opt$tissue_index]
}

message("Running tissue: ", tissue_to_run)
df_sub <- df[df$tissue == tissue_to_run, ]

loci <- find_loci_twas(df_sub)
loci$pip[is.na(loci$pip)] <- 0
loci <- loci[loci$pip >= opt$pip, ]
if (nrow(loci) == 0) {
  message("No TWAS genes passed PIP threshold. Nothing to write.")
  quit(save = "no", status = 0)
}

chunks <- parse_chunks(opt$chunks)
message("Reading processed CMap chunks: ", paste(chunks, collapse = ", "))

cmap_data <- data.frame()
for (chunk in chunks) {
  file_i <- fill_template(opt$processed_cmap_template, opt$cell, chunk)
  check_file(file_i, paste0("processed CMap chunk ", chunk))
  tmp <- fread(file_i)
  if (ncol(tmp) < 2) stop("Processed CMap file has fewer than 2 columns: ", file_i, call. = FALSE)

  idx <- match(loci$entrez_id, tmp[[1]])
  tmp_sub <- tmp[idx, -1, with = FALSE]
  chunk_mat <- as.data.frame(t(tmp_sub))
  rownames(chunk_mat) <- names(tmp)[-1]
  cmap_data <- rbind(cmap_data, chunk_mat)
}

colnames(cmap_data) <- loci$entrez_id
if (ncol(cmap_data) <= 1) {
  message("Fewer than two TWAS genes matched processed CMap data. Nothing to write.")
  quit(save = "no", status = 0)
}

cmap_data <- as.data.frame(lapply(cmap_data, as.numeric), row.names = rownames(cmap_data))
cmap_data[abs(cmap_data) <= opt$z_threshold] <- 0
row_signal <- rowSums(abs(cmap_data), na.rm = TRUE)
cmap_data <- cmap_data[row_signal > 0, , drop = FALSE]

pert_id <- vapply(strsplit(rownames(cmap_data), ":"), function(x) x[1], character(1))
cmap_data <- cmap_data[pert_id %in% drug_sub$pert_id, , drop = FALSE]
if (nrow(cmap_data) == 0) {
  message("No processed CMap perturbations matched validation drugs. Nothing to write.")
  quit(save = "no", status = 0)
}

message("Running bootstrap with ", opt$n_bootstrap, " replicates...")
set.seed(opt$seed)

data_bootstrap <- vector("list", ncol(cmap_data))
for (j in seq_len(ncol(cmap_data))) {
  set.seed(opt$seed + j)
  signal_j <- matrix(cmap_data[[j]], nrow = nrow(cmap_data), ncol = opt$n_bootstrap)
  noise_j <- matrix(rnorm(opt$n_bootstrap * nrow(cmap_data), mean = 0, sd = sd(cmap_data[[j]], na.rm = TRUE)),
                    nrow = nrow(cmap_data), ncol = opt$n_bootstrap)
  data_bootstrap[[j]] <- signal_j + noise_j
}

res <- sapply(seq_len(opt$n_bootstrap), function(x) {
  data_temp <- cmap_data
  for (j in seq_len(ncol(cmap_data))) {
    data_temp[[j]] <- data_bootstrap[[j]][, x]
  }
  data_sd <- sqrt(rowSums(as.matrix(data_temp)^2, na.rm = TRUE))
  data_z <- as.matrix(data_temp) %*% as.matrix(loci$TWAS.Z)
  as.numeric(data_z / data_sd)
})

data_sd <- sqrt(rowSums(as.matrix(cmap_data)^2, na.rm = TRUE))
data_z <- as.matrix(cmap_data) %*% as.matrix(loci$TWAS.Z)

out <- data.frame(
  drug = rownames(cmap_data),
  tissue = tissue_to_run,
  zscore = as.numeric(data_z / data_sd),
  zmean = rowMeans(res, na.rm = TRUE),
  zsd = apply(res, 1, sd, na.rm = TRUE),
  stringsAsFactors = FALSE
)
out$pval <- pchisq((out$zscore / out$zsd)^2, df = 1, lower.tail = FALSE)

safe_tissue <- gsub("[^A-Za-z0-9_.-]+", "_", tissue_to_run)
safe_cell <- gsub("[^A-Za-z0-9_.-]+", "_", opt$cell)
safe_model <- gsub("[^A-Za-z0-9_.-]+", "_", opt$model_name)
outfile <- file.path(
  opt$output_dir,
  paste0(opt$output_prefix, ".", safe_model, "_", safe_tissue, "_", safe_cell, "_PIP_", opt$pip, ".txt")
)

fwrite(out, outfile, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
message("Wrote: ", outfile)
