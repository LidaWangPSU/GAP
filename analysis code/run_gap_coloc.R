#!/usr/bin/env Rscript

# Colocalization analysis for GAP loci and external molecular QTL datasets
#
# GitHub/shared version:
#   - No private hard-coded paths
#   - All input/output locations are supplied as command-line arguments
#   - Dataset selection uses source names rather than fragile numeric ranges where possible
#
# Example:
# Rscript run_gap_coloc.R \
#   --eqtl-index=1 \
#   --data-name=BIOVU \
#   --gap-model-dir=/path/to/gap_model \
#   --loci-rds=/path/to/BIOVU_union+100snps-5e-04-lrt_gc.rds \
#   --gtex-eqtl-dir=/path/to/GTEx_Analysis_v7_eQTL_all_associations \
#   --geuvadis-dir=/path/to/GEUVADIS/Cells_EBV-transformed_lymphocytes \
#   --eqtlgen-dir=/path/to/eQTLGenA_new \
#   --jobs-eqtl-dir=/path/to/JOBS/eQTLGen_1k1k/eqtls \
#   --sle-eqtl-dir=/path/to/sle_eqtls/eQTLs \
#   --gencode-gtf=/path/to/gencode.v34.GRCh38.genes.collapsed_only.gtf \
#   --caqtl-dir=/path/to/caqtl/Kumasaka \
#   --mqtl-dir=/path/to/GTEx_mqtls/mQTLs/WholeBlood_bychr \
#   --gtex-sqtl-dir=/path/to/GTEx_sQTL/GTEx_v8/sQTL_sumstat \
#   --aida-nominal-dir=/path/to/AIDA/nominal_pass \
#   --liftover-chain=/path/to/hg38ToHg19.over.chain \
#   --mixqtl-dir=/path/to/GTEx_mixQTL \
#   --pqtl-dir=/path/to/pqtl \
#   --sle-eqtl-discov-dir=/path/to/sle_eqtls/eQTLs_dis_cov \
#   --output-dir=/path/to/output

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
  library(rtracklayer)
  library(GenomicRanges)
})

parse_args <- function() {
  raw <- commandArgs(trailingOnly = TRUE)

  # Backward-compatible mode: first positional argument is eqtl index.
  if (length(raw) == 1 && !grepl("^--", raw[1])) {
    return(list(eqtl_index = as.integer(raw[1])))
  }

  args <- list()
  for (x in raw) {
    if (!grepl("^--", x)) next
    kv <- strsplit(sub("^--", "", x), "=", fixed = TRUE)[[1]]
    key <- gsub("-", "_", kv[1])
    value <- if (length(kv) > 1) paste(kv[-1], collapse = "=") else TRUE
    args[[key]] <- value
  }

  if (!is.null(args$eqtl_index)) args$eqtl_index <- as.integer(args$eqtl_index)
  args
}

get_arg <- function(args, name, default = NULL, required = FALSE) {
  value <- args[[name]]
  if (is.null(value) || is.na(value) || identical(value, "")) {
    if (required) stop("Missing required argument: --", gsub("_", "-", name), call. = FALSE)
    return(default)
  }
  value
}

assert_file <- function(path, label) {
  if (is.null(path) || !file.exists(path)) stop(label, " does not exist: ", path, call. = FALSE)
  invisible(path)
}

assert_dir <- function(path, label, required = TRUE) {
  if (is.null(path) || !dir.exists(path)) {
    if (required) stop(label, " does not exist: ", path, call. = FALSE)
    return(NULL)
  }
  invisible(path)
}

read_required <- function(path, ...) {
  assert_file(path, "Input file")
  fread(path, ...)
}

make_variant_id <- function(chr, pos, ref, alt, build = "b37") {
  paste(chr, pos, ref, alt, build, sep = "_")
}

safe_min_nonzero <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  nz <- x[x != 0]
  if (length(nz) == 0) return(min(x, na.rm = TRUE))
  min(nz, na.rm = TRUE)
}

load_gencode_gene_map <- function(gtf_path) {
  assert_file(gtf_path, "Gencode GTF")
  hg38_gtf <- rtracklayer::import(gtf_path)
  hg38_gtf <- hg38_gtf[hg38_gtf$type == "gene"]
  data.frame(
    chr = as.character(seqnames(hg38_gtf)),
    start = start(hg38_gtf),
    end = end(hg38_gtf),
    gene_length = width(hg38_gtf),
    gene_id_full = hg38_gtf$gene_id,
    gene_type = hg38_gtf$gene_type,
    gene_name = hg38_gtf$gene_name,
    stringsAsFactors = FALSE
  )
}

load_eqtl <- function(source, chr, gwas_sub, cfg, source_sets, gencode_map = NULL, chain = NULL) {
  tmp <- NULL
  N_eqtl <- NA_real_

  if (source == "GEUVADIS") {
    path <- file.path(cfg$geuvadis_dir, paste0("GEUVADIS_chr", chr, "_Cells_EBV-transformed_lymphocytes_sumstat.txt"))
    tmp <- as.data.frame(read_required(path))
    snps <- intersect(gwas_sub$snp, tmp$V11)
    tmp <- tmp[tmp$V11 %in% snps, , drop = FALSE]
    tmp$pval <- pchisq((tmp$V8 / tmp$V9)^2, df = 1, lower.tail = FALSE)
    tmp <- tmp[, c("V1", "V11", "pval")]
    colnames(tmp) <- c("gene", "snp", "pval")
    N_eqtl <- 358

  } else if (source == "eQTLGen") {
    path <- file.path(cfg$eqtlgen_dir, paste0("eQTLGenA_Whole_Blood_chr", chr, ".txt"))
    tmp <- as.data.frame(read_required(path))
    tmp$rs_flip <- make_variant_id(tmp$chr, tmp$pos, tmp$alt, tmp$ref)
    snps <- intersect(gwas_sub$snp, tmp$rs)
    snps2 <- intersect(gwas_sub$snp, tmp$rs_flip)
    tmp <- tmp[tmp$rs %in% snps | tmp$rs_flip %in% snps2, , drop = FALSE]
    tmp$pval <- pchisq((tmp$beta.est / tmp$beta.se)^2, df = 1, lower.tail = FALSE)
    tmp1 <- tmp[tmp$rs %in% snps, c("gene", "rs", "pval"), drop = FALSE]
    tmp2 <- tmp[tmp$rs_flip %in% snps2, c("gene", "rs_flip", "pval"), drop = FALSE]
    colnames(tmp2) <- colnames(tmp1)
    tmp <- rbind(tmp1, tmp2)
    colnames(tmp) <- c("gene", "snp", "pval")
    N_eqtl <- 25000

  } else if (source %in% source_sets$onek1k_cell_types) {
    path <- file.path(cfg$jobs_eqtl_dir, paste0(source, "_JOBS_summary_statistics.txt"))
    tmp <- as.data.frame(read_required(path))
    snps <- intersect(gwas_sub$snp, tmp$snp_hg19)
    tmp <- tmp[tmp$snp_hg19 %in% snps, , drop = FALSE]
    tmp$pval <- pchisq((tmp$beta.est / tmp$beta.se)^2, df = 1, lower.tail = FALSE)
    tmp <- tmp[, c("gene", "snp_hg19", "pval")]
    colnames(tmp) <- c("gene", "snp", "pval")
    N_eqtl <- 1000

  } else if (source %in% source_sets$gtex_tissues) {
    path <- file.path(cfg$gtex_eqtl_dir, paste0(source, ".allpairs.txt.gz"))
    tmp <- as.data.frame(read_required(path))
    snps <- intersect(gwas_sub$snp, tmp$variant_id)
    tmp <- tmp[tmp$variant_id %in% snps, , drop = FALSE]
    tmp$pval <- pchisq((tmp$slope / tmp$slope_se)^2, df = 1, lower.tail = FALSE)
    tmp$gene_id <- sub("\\..*$", "", tmp$gene_id)
    N_eqtl <- median(tmp$ma_samples, na.rm = TRUE)
    tmp <- tmp[, c("gene_id", "variant_id", "pval")]
    colnames(tmp) <- c("gene", "snp", "pval")

  } else if (source %in% source_sets$sle_cells) {
    if (is.null(gencode_map)) stop("--gencode-gtf is required for SLE eQTL sources", call. = FALSE)
    path <- file.path(cfg$sle_eqtl_dir, paste0(source, "_cis_eQTLs.txt"))
    tmp <- as.data.frame(read_required(path))
    tmp$gene_id <- gencode_map$gene_id_full[match(tmp$gene, gencode_map$gene_name)]
    tmp$gene_id <- sub("\\..*$", "", tmp$gene_id)
    tmp$rs_flip <- make_variant_id(tmp$chr, tmp$pos, tmp$alt, tmp$ref)
    snps <- intersect(gwas_sub$snp, tmp$rs)
    snps2 <- intersect(gwas_sub$snp, tmp$rs_flip)
    tmp <- tmp[tmp$rs %in% snps | tmp$rs_flip %in% snps2, , drop = FALSE]
    tmp$pval <- pchisq((tmp$beta.est / tmp$beta.se)^2, df = 1, lower.tail = FALSE)
    N_eqtl <- ceiling(median(tmp$n, na.rm = TRUE))
    tmp1 <- tmp[tmp$rs %in% snps, c("gene_id", "rs", "pval"), drop = FALSE]
    tmp2 <- tmp[tmp$rs_flip %in% snps2, c("gene_id", "rs_flip", "pval"), drop = FALSE]
    colnames(tmp2) <- colnames(tmp1)
    tmp <- rbind(tmp1, tmp2)
    colnames(tmp) <- c("gene", "snp", "pval")

  } else if (source == "caQTL") {
    path <- file.path(cfg$caqtl_dir, paste0("all.chr", chr, ".caQTL.summary.statistics"))
    tmp <- as.data.frame(read_required(path))
    tmp$gene <- paste0("Peak_", tmp$V1)
    tmp$rs <- make_variant_id(tmp$V2, tmp$V3, tmp$V5, tmp$V6)
    snps <- intersect(gwas_sub$snp, tmp$rs)
    tmp <- tmp[tmp$rs %in% snps, , drop = FALSE]
    tmp$pval <- pchisq((tmp$V10 / tmp$V11)^2, df = 1, lower.tail = FALSE)
    tmp <- tmp[, c("gene", "rs", "pval")]
    colnames(tmp) <- c("gene", "snp", "pval")
    N_eqtl <- 100

  } else if (source == "mQTL") {
    path <- file.path(cfg$mqtl_dir, paste0("chr", chr, "_liftover.txt"))
    tmp <- as.data.frame(read_required(path))
    tmp <- tmp[tmp$V3 %in% gwas_sub$pos, , drop = FALSE]
    tmp$rs <- make_variant_id(gsub("chr", "", tmp$V1), tmp$V3, tmp$V15, tmp$V16)
    snps <- intersect(gwas_sub$snp, tmp$rs)
    tmp <- tmp[tmp$rs %in% snps, c("V6", "rs", "V12"), drop = FALSE]
    colnames(tmp) <- c("gene", "snp", "pval")
    N_eqtl <- 500

  } else if (source %in% source_sets$gtex_sqtl_sources) {
    tissue <- sub("^sQTL_", "", source)
    path <- file.path(cfg$gtex_sqtl_dir, tissue, paste0(tissue, "_chr", chr, "_ALL_cis_sQTLs.txt"))
    tmp <- as.data.frame(read_required(path))
    tmp <- tmp[tmp$V10 %in% gwas_sub$pos, , drop = FALSE]
    tmp$rs <- make_variant_id(tmp$V2, tmp$V10, tmp$V4, tmp$V5)
    tmp$rs_flip <- make_variant_id(tmp$V2, tmp$V10, tmp$V5, tmp$V4)
    snps <- intersect(gwas_sub$snp, tmp$rs)
    snps2 <- intersect(gwas_sub$snp, tmp$rs_flip)
    tmp <- tmp[tmp$rs %in% snps | tmp$rs_flip %in% snps2, , drop = FALSE]
    tmp$pval <- pchisq((tmp$V6 / tmp$V7)^2, df = 1, lower.tail = FALSE)
    N_eqtl <- ceiling(median(tmp$V8, na.rm = TRUE))
    tmp1 <- tmp[tmp$rs %in% snps, c("V1", "rs", "pval"), drop = FALSE]
    tmp2 <- tmp[tmp$rs_flip %in% snps2, c("V1", "rs_flip", "pval"), drop = FALSE]
    colnames(tmp2) <- colnames(tmp1)
    tmp <- rbind(tmp1, tmp2)
    colnames(tmp) <- c("gene", "snp", "pval")

  } else if (source %in% source_sets$aida_sqtl_sources) {
    if (is.null(chain)) stop("--liftover-chain is required for AIDA sQTL sources", call. = FALSE)
    tissue <- sub("^sQTL_", "", source)
    path <- file.path(cfg$aida_nominal_dir, tissue, paste0(tissue, "_", chr, "_nominal_1.txt.gz"))
    tmp <- as.data.frame(read_required(path, sep = " ", header = FALSE))
    tmp <- tmp[tmp$V12 == tmp$V13, , drop = FALSE]
    if (nrow(tmp) == 0) return(NULL)

    tmp$start <- tmp$V13 - 1
    tmp$end <- tmp$V13
    tmp$chr <- tmp$V11
    gr <- makeGRangesFromDataFrame(tmp, keep.extra.columns = TRUE)
    tmp <- as.data.frame(unlist(liftOver(gr, chain)), row.names = NULL)
    tmp <- tmp[tmp$end %in% gwas_sub$pos, , drop = FALSE]
    if (nrow(tmp) == 0) return(NULL)

    tmp$ref <- vapply(strsplit(tmp$V10, ":"), `[`, character(1), 3)
    tmp$alt <- vapply(strsplit(tmp$V10, ":"), `[`, character(1), 4)
    tmp$seqnames <- gsub("chr", "", tmp$seqnames)
    tmp$rs <- make_variant_id(tmp$seqnames, tmp$end, tmp$ref, tmp$alt)
    tmp$rs_flip <- make_variant_id(tmp$seqnames, tmp$end, tmp$alt, tmp$ref)
    snps <- intersect(gwas_sub$snp, tmp$rs)
    snps2 <- intersect(gwas_sub$snp, tmp$rs_flip)
    tmp <- tmp[tmp$rs %in% snps | tmp$rs_flip %in% snps2, , drop = FALSE]
    tmp$pair <- paste0(tmp$V6, tmp$V1)
    N_eqtl <- cfg$aida_n_eqtl
    tmp1 <- tmp[tmp$rs %in% snps, c("pair", "rs", "V14"), drop = FALSE]
    tmp2 <- tmp[tmp$rs_flip %in% snps2, c("pair", "rs_flip", "V14"), drop = FALSE]
    colnames(tmp2) <- colnames(tmp1)
    tmp <- rbind(tmp1, tmp2)
    colnames(tmp) <- c("gene", "snp", "pval")

  } else if (source %in% source_sets$mixqtl_sources) {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      stop("Package 'arrow' is required for mixQTL sources", call. = FALSE)
    }
    if (is.null(chain)) stop("--liftover-chain is required for mixQTL sources", call. = FALSE)
    tissue <- sub("^mixQTL_", "", source)
    path <- file.path(cfg$mixqtl_dir, tissue, paste0("mixqtl.", tissue, ".chr", chr, ".parquet"))
    assert_file(path, "mixQTL parquet file")
    tmp <- as.data.frame(arrow::read_parquet(path))
    tmp <- na.omit(tmp[, c("phenotype_id", "variant_id", "pval_asc")])
    tmp$chr <- vapply(strsplit(tmp$variant_id, "_"), `[`, character(1), 1)
    tmp$end <- as.numeric(vapply(strsplit(tmp$variant_id, "_"), `[`, character(1), 2))
    tmp$start <- tmp$end - 1

    gr <- makeGRangesFromDataFrame(tmp, keep.extra.columns = TRUE)
    tmp <- as.data.frame(unlist(liftOver(gr, chain)), row.names = NULL)
    tmp <- tmp[tmp$end %in% gwas_sub$pos, , drop = FALSE]
    if (nrow(tmp) == 0) return(NULL)

    tmp$ref <- vapply(strsplit(tmp$variant_id, "_"), `[`, character(1), 3)
    tmp$alt <- vapply(strsplit(tmp$variant_id, "_"), `[`, character(1), 4)
    tmp$seqnames <- gsub("chr", "", tmp$seqnames)
    tmp$rs <- make_variant_id(tmp$seqnames, tmp$end, tmp$ref, tmp$alt)
    tmp$rs_flip <- make_variant_id(tmp$seqnames, tmp$end, tmp$alt, tmp$ref)
    snps <- intersect(gwas_sub$snp, tmp$rs)
    snps2 <- intersect(gwas_sub$snp, tmp$rs_flip)
    tmp <- tmp[tmp$rs %in% snps | tmp$rs_flip %in% snps2, , drop = FALSE]
    N_eqtl <- cfg$mixqtl_n_eqtl
    tmp1 <- tmp[tmp$rs %in% snps, c("phenotype_id", "rs", "pval_asc"), drop = FALSE]
    tmp2 <- tmp[tmp$rs_flip %in% snps2, c("phenotype_id", "rs_flip", "pval_asc"), drop = FALSE]
    colnames(tmp2) <- colnames(tmp1)
    tmp <- rbind(tmp1, tmp2)
    colnames(tmp) <- c("gene", "snp", "pval")

  } else if (source %in% source_sets$pqtl_sources) {
    path <- file.path(cfg$pqtl_dir, paste0("processed_", source, ".txt"))
    tmp <- read_required(path)
    tmp <- tmp[, 1:3]
    snps <- intersect(gwas_sub$snp, tmp$snp)
    tmp <- as.data.frame(tmp[tmp$snp %in% snps, ])
    colnames(tmp)[1:3] <- c("gene", "snp", "pval")
    N_eqtl <- cfg$pqtl_n_eqtl

  } else if (source %in% source_sets$sle_cells_new) {
    source_file <- gsub("NEW", "ALL", source)
    path <- file.path(cfg$sle_eqtl_discov_dir, paste0(source_file, "_chr", chr, "_cis_eQTLs.txt"))
    tmp <- as.data.frame(read_required(path))
    tmp$gene <- sub("\\..*$", "", tmp$gene)
    tmp$rs_flip <- make_variant_id(tmp$chr, tmp$pos, tmp$alt, tmp$ref)
    snps <- intersect(gwas_sub$snp, tmp$rs)
    snps2 <- intersect(gwas_sub$snp, tmp$rs_flip)
    tmp <- tmp[tmp$rs %in% snps | tmp$rs_flip %in% snps2, , drop = FALSE]
    tmp$pval <- pchisq((tmp$beta.est / tmp$beta.se)^2, df = 1, lower.tail = FALSE)
    N_eqtl <- ceiling(median(tmp$n, na.rm = TRUE))
    tmp1 <- tmp[tmp$rs %in% snps, c("gene", "rs", "pval"), drop = FALSE]
    tmp2 <- tmp[tmp$rs_flip %in% snps2, c("gene", "rs_flip", "pval"), drop = FALSE]
    colnames(tmp2) <- colnames(tmp1)
    tmp <- rbind(tmp1, tmp2)
    colnames(tmp) <- c("gene", "snp", "pval")

  } else {
    stop("Unsupported QTL source: ", source, call. = FALSE)
  }

  if (is.null(tmp) || nrow(tmp) == 0) return(NULL)
  tmp <- tmp[!is.na(tmp$gene) & !is.na(tmp$snp) & !is.na(tmp$pval), , drop = FALSE]
  if (nrow(tmp) == 0) return(NULL)

  list(data = tmp, n_eqtl = N_eqtl)
}

main <- function() {
  args <- parse_args()

  cfg <- list(
    data_name = get_arg(args, "data_name", "BIOVU"),
    eqtl_index = get_arg(args, "eqtl_index", required = TRUE),

    gap_model_dir = get_arg(args, "gap_model_dir", required = TRUE),
    loci_rds = get_arg(args, "loci_rds", required = TRUE),
    output_dir = get_arg(args, "output_dir", required = TRUE),

    gtex_eqtl_dir = get_arg(args, "gtex_eqtl_dir", NULL),
    geuvadis_dir = get_arg(args, "geuvadis_dir", NULL),
    eqtlgen_dir = get_arg(args, "eqtlgen_dir", NULL),
    jobs_eqtl_dir = get_arg(args, "jobs_eqtl_dir", NULL),
    sle_eqtl_dir = get_arg(args, "sle_eqtl_dir", NULL),
    gencode_gtf = get_arg(args, "gencode_gtf", NULL),
    caqtl_dir = get_arg(args, "caqtl_dir", NULL),
    mqtl_dir = get_arg(args, "mqtl_dir", NULL),
    gtex_sqtl_dir = get_arg(args, "gtex_sqtl_dir", NULL),
    aida_nominal_dir = get_arg(args, "aida_nominal_dir", NULL),
    liftover_chain = get_arg(args, "liftover_chain", NULL),
    mixqtl_dir = get_arg(args, "mixqtl_dir", NULL),
    pqtl_dir = get_arg(args, "pqtl_dir", NULL),
    sle_eqtl_discov_dir = get_arg(args, "sle_eqtl_discov_dir", NULL),

    aida_n_eqtl = as.numeric(get_arg(args, "aida_n_eqtl", 474)),
    mixqtl_n_eqtl = as.numeric(get_arg(args, "mixqtl_n_eqtl", 500)),
    pqtl_n_eqtl = as.numeric(get_arg(args, "pqtl_n_eqtl", 50000)),
    min_eqtl_snps = as.integer(get_arg(args, "min_eqtl_snps", 100)),
    coloc_window_bp = as.integer(get_arg(args, "coloc_window_bp", 1000000)),
    gwas_p_threshold = as.numeric(get_arg(args, "gwas_p_threshold", 5e-8))
  )

  dir.create(cfg$output_dir, showWarnings = FALSE, recursive = TRUE)

  source_sets <- list()
  source_sets$onek1k_cell_types <- c("B_IN", "B_Mem", "Plasma", "CD4_ET", "CD4_NC", "CD4_SOX4",
                                     "CD8_NC", "CD8_ET", "CD8_S100B", "DC1", "Mono_NC",
                                     "Mono_C", "NK", "NK_R")

  if (!is.null(cfg$gtex_eqtl_dir) && dir.exists(cfg$gtex_eqtl_dir)) {
    source_sets$gtex_tissues <- sub("\\.allpairs\\.txt\\.gz$", "", basename(dir(cfg$gtex_eqtl_dir, pattern = "\\.allpairs\\.txt\\.gz$")))
  } else {
    source_sets$gtex_tissues <- character(0)
  }

  source_sets$sle_cells <- paste0(c("B", "CD8", "Progen", "cM", "NK", "Prolif", "ncM", "CD4", "PB", "cDC", "pDC"), "_ALL")
  source_sets$sle_cells_new <- paste0(c("B", "T8", "Progen", "cM", "NK", "Prolif", "ncM", "T4", "PB", "cDC", "pDC"), "_NEW")

  source_sets$gtex_sqtl_sources <- c("sQTL_Whole_Blood", "sQTL_Cells_EBV-transformed_lymphocytes")

  if (!is.null(cfg$aida_nominal_dir) && dir.exists(cfg$aida_nominal_dir)) {
    source_sets$aida_sqtl_sources <- paste0("sQTL_", basename(dir(cfg$aida_nominal_dir)))
  } else {
    source_sets$aida_sqtl_sources <- character(0)
  }

  source_sets$mixqtl_sources <- c("mixQTL_Whole_Blood", "mixQTL_Cells_EBV-transformed_lymphocytes")
  source_sets$pqtl_sources <- c(paste0("olink_100_", 1:100), paste0("spmascan_100_", 1:100))

  eqtls_data <- c(
    "GEUVADIS",
    "eQTLGen",
    source_sets$onek1k_cell_types,
    source_sets$gtex_tissues,
    source_sets$sle_cells,
    "caQTL",
    "mQTL",
    source_sets$gtex_sqtl_sources,
    source_sets$aida_sqtl_sources,
    source_sets$mixqtl_sources,
    source_sets$pqtl_sources,
    source_sets$sle_cells_new
  )

  if (cfg$eqtl_index < 1 || cfg$eqtl_index > length(eqtls_data)) {
    stop("--eqtl-index is out of range. Valid range: 1-", length(eqtls_data), call. = FALSE)
  }
  qtl_source <- eqtls_data[cfg$eqtl_index]
  message("Selected QTL source: ", qtl_source)

  gencode_map <- NULL
  if (qtl_source %in% source_sets$sle_cells) {
    gencode_map <- load_gencode_gene_map(cfg$gencode_gtf)
  }

  chain <- NULL
  if (qtl_source %in% c(source_sets$aida_sqtl_sources, source_sets$mixqtl_sources)) {
    assert_file(cfg$liftover_chain, "LiftOver chain")
    chain <- import.chain(cfg$liftover_chain)
  }

  # Load GAP model output across chromosomes.
  gwas_list <- vector("list", 22)
  for (chr in 1:22) {
    path <- file.path(
      cfg$gap_model_dir,
      paste0("GAP_", cfg$data_name, "_union+100snps-5e-04_chr", chr, "_lrt_new.txt")
    )
    gwas_list[[chr]] <- read_required(path)
  }
  df <- as.data.frame(rbindlist(gwas_list, fill = TRUE))
  variant_parts <- tstrsplit(df$snp, "_")
  df$chr <- as.numeric(variant_parts[[1]])
  df$pos <- as.numeric(variant_parts[[2]])
  df$ref <- as.character(variant_parts[[3]])
  df$alt <- as.character(variant_parts[[4]])

  loci <- readRDS(assert_file(cfg$loci_rds, "Loci RDS"))
  loci01 <- loci[[1]]
  loci12 <- loci[[2]]
  loci01$type <- "01"
  loci12$type <- "12"
  colnames(loci12) <- colnames(loci01)
  loci <- rbind(loci01, loci12)
  loci$pos <- as.numeric(vapply(strsplit(loci$snp, "_"), `[`, character(1), 2))

  out <- list()
  out_i <- 0L

  for (l in seq_len(nrow(loci))) {
    chr <- loci$chr[l]
    pos <- loci$pos[l]
    type <- loci$type[l]

    gwas_sub <- df[df$chr == chr & df$pos < (pos + cfg$coloc_window_bp) & df$pos > (pos - cfg$coloc_window_bp), , drop = FALSE]
    if (nrow(gwas_sub) == 0) next

    N_gwas <- if (type == "01") 10503 else if (type == "12") 2361 else NA_real_

    eqtl_obj <- tryCatch(
      load_eqtl(qtl_source, chr, gwas_sub, cfg, source_sets, gencode_map = gencode_map, chain = chain),
      error = function(e) {
        message("Skipping locus ", loci$snp[l], " because QTL loading failed: ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(eqtl_obj)) next

    tmp <- eqtl_obj$data
    N_eqtl <- eqtl_obj$n_eqtl
    genes <- unique(tmp$gene)

    for (g in genes) {
      eqtls <- tmp[tmp$gene == g, , drop = FALSE]
      if (nrow(eqtls) < cfg$min_eqtl_snps) next

      input <- merge(gwas_sub, eqtls, by = "snp")
      if (nrow(input) == 0) next

      if (type == "01") {
        input$pval_gwas <- pchisq(input$zscore_01^2, df = 1, lower.tail = FALSE)
      } else if (type == "12") {
        input$pval_gwas <- pchisq(input$zscore_12^2, df = 1, lower.tail = FALSE)
      } else {
        next
      }

      if (all(is.na(input$pval_gwas)) || min(input$pval_gwas, na.rm = TRUE) > cfg$gwas_p_threshold) next

      input$pval[input$pval == 0] <- safe_min_nonzero(input$pval)
      input$pval_gwas[input$pval_gwas == 0] <- safe_min_nonzero(input$pval_gwas)
      input$MAF <- 0.5

      dataset1 <- list(snp = input$snp, pvalues = input$pval_gwas, type = "quant", N = N_gwas)
      dataset2 <- list(snp = input$snp, pvalues = input$pval, type = "quant", sdY = 1, N = N_eqtl)

      coloc_fit <- tryCatch(
        coloc.abf(dataset1, dataset2, MAF = input$MAF),
        error = function(e) NULL
      )
      if (is.null(coloc_fit)) next

      out_i <- out_i + 1L
      out[[out_i]] <- data.frame(
        data = qtl_source,
        loci = loci$snp[l],
        type = type,
        chr = chr,
        matched_snp = nrow(input),
        PP4 = as.numeric(coloc_fit$summary["PP.H4.abf"]),
        gene = g,
        min_eqtls_p = min(input$pval, na.rm = TRUE),
        min_gwas_p = min(input$pval_gwas, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }

  out_df <- if (length(out) > 0) rbindlist(out, fill = TRUE) else data.table()
  message("Number of colocalization results: ", nrow(out_df))
  if (nrow(out_df) > 0) print(summary(as.numeric(out_df$PP4)))

  output_file <- file.path(cfg$output_dir, paste0(cfg$data_name, "_GAP_", qtl_source, "_coloc.txt"))
  fwrite(out_df, output_file, sep = "\t", quote = FALSE, na = "NA")
  message("Saved: ", output_file)
}

main()
