#!/usr/bin/env Rscript

# Run GAP TWAS with external expression-prediction weights and PLINK LD reference.
#
# This GitHub-shareable version removes user-specific/private paths from the
# original analysis script. All input locations are supplied as command-line
# arguments.
#
# Required R packages:
#   install.packages(c("data.table", "optparse", "BEDMatrix"))
#
# Example:
# Rscript run_gap_twas_github_shared.R \
#   --model union+100snps-5e-04-lrt \
#   --stage 01 \
#   --data-name BIOVU \
#   --weights-file /path/to/EXPRESSO_eQTLGen_weights_hg19.txt \
#   --weights-name eQTLGen \
#   --ldref-prefix-template /path/to/1000G_hg19_0.05/EUR.chr{chr}.final \
#   --gap-model-dir /path/to/gap_model \
#   --output-dir /path/to/twas_output
#
# Model-specific inputs:
#   model = mtag:
#     --mtag-dir must contain Sumstat_<data-name>_trait_1.txt and
#     Sumstat_<data-name>_trait_2.txt.
#   model = gap-trans:
#     provide --gap-trans-stage01-file and/or --gap-trans-stage12-file.
#   model = no-prior-lrt:
#     --mixture-model-dir must contain GAP_mixture_<data-name>_chr<chr>_lrt.txt.
#   model = union+100snps-5e-04-lrt:
#     --gap-model-dir must contain GAP_<data-name>_union+100snps-5e-04_chr<chr>_lrt_new.txt.

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(BEDMatrix)
})

options(stringsAsFactors = FALSE)

option_list <- list(
  make_option("--model", type = "character", default = "union+100snps-5e-04-lrt",
              help = "GWAS model: mtag, no-prior-lrt, union+100snps-5e-04-lrt, or gap-trans [default: %default]"),
  make_option("--stage", type = "character", default = "01",
              help = "Disease-progression stage: 01 or 12 [default: %default]"),
  make_option("--data-name", type = "character", default = "BIOVU",
              help = "Dataset prefix used in GAP filenames [default: %default]"),
  make_option("--weights-file", type = "character", default = NULL,
              help = "Expression-prediction weight file. Required."),
  make_option("--weights-name", type = "character", default = NULL,
              help = "Short name for the weight file, used in output filename. If omitted, basename(weights-file) is used."),
  make_option("--ldref-prefix-template", type = "character", default = NULL,
              help = "PLINK LD reference prefix template with {chr}, e.g. /path/EUR.chr{chr}.final. Required."),
  make_option("--mtag-dir", type = "character", default = NULL,
              help = "Directory containing MTAG summary statistics."),
  make_option("--gap-trans-stage01-file", type = "character", default = NULL,
              help = "Stage 01 trans/meta summary statistics file for model=gap-trans."),
  make_option("--gap-trans-stage12-file", type = "character", default = NULL,
              help = "Stage 12 trans/meta summary statistics file for model=gap-trans."),
  make_option("--mixture-model-dir", type = "character", default = NULL,
              help = "Directory containing GAP_mixture_<data-name>_chr<chr>_lrt.txt files."),
  make_option("--gap-model-dir", type = "character", default = NULL,
              help = "Directory containing GAP_<data-name>_union+100snps-5e-04_chr<chr>_lrt_new.txt files."),
  make_option("--output-dir", type = "character", default = ".",
              help = "Output directory [default: %default]"),
  make_option("--min-r2pred", type = "double", default = 0.6,
              help = "Minimum mean imputation r2 at non-zero expression-weight SNPs [default: %default]"),
  make_option("--max-impute", type = "double", default = 0.5,
              help = "Maximum fraction of missing GWAS Z-scores allowed per gene [default: %default]"),
  make_option("--ld-ridge", type = "double", default = 0.1,
              help = "Small ridge added to the LD diagonal for numerical stability [default: %default]"),
  make_option("--chromosomes", type = "character", default = "1:22",
              help = "Chromosomes to analyze, e.g. 1:22 or 1,2,3 [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

parse_chromosomes <- function(x) {
  if (grepl(":", x, fixed = TRUE)) {
    z <- as.integer(strsplit(x, ":", fixed = TRUE)[[1]])
    return(seq(z[1], z[2]))
  }
  as.integer(strsplit(x, ",", fixed = TRUE)[[1]])
}

stop_if_missing <- function(path, label = "file") {
  if (is.null(path) || is.na(path) || !file.exists(path)) {
    stop(sprintf("Missing %s: %s", label, path), call. = FALSE)
  }
}

stop_if_missing_dir <- function(path, label = "directory") {
  if (is.null(path) || is.na(path) || !dir.exists(path)) {
    stop(sprintf("Missing %s: %s", label, path), call. = FALSE)
  }
}

plink_prefix_for_chr <- function(template, chr) {
  gsub("\\{chr\\}", chr, template)
}

read_plink_custom <- function(root, impute = c("none", "avg")) {
  impute <- match.arg(impute)
  bedfile <- paste0(path.expand(root), ".bed")
  famfile <- paste0(path.expand(root), ".fam")
  bimfile <- paste0(path.expand(root), ".bim")
  stop_if_missing(bedfile, "PLINK .bed file")
  stop_if_missing(famfile, "PLINK .fam file")
  stop_if_missing(bimfile, "PLINK .bim file")

  bim <- read.table(bimfile, header = FALSE, sep = "", stringsAsFactors = FALSE)
  fam <- read.table(famfile, header = FALSE, sep = "", stringsAsFactors = FALSE)
  geno <- BEDMatrix::BEDMatrix(bedfile, n = nrow(fam), p = nrow(bim))
  geno <- as.matrix(geno)

  if (impute == "avg") {
    geno_na <- is.na(geno)
    if (any(geno_na)) {
      means <- colMeans(geno, na.rm = TRUE)
      geno[geno_na] <- rep(means, colSums(geno_na))
    }
  }

  colnames(geno) <- bim[, 2]
  rownames(geno) <- paste(fam[, 1], fam[, 2], sep = ":")
  list(bed = geno, fam = fam, bim = bim)
}

harmonize_ld_reference <- function(genos) {
  keep <- !duplicated(genos$bim[, 4])
  genos$bim <- genos$bim[keep, ]
  genos$bed <- genos$bed[, genos$bim[, 2], drop = FALSE]
  genos$bim$V2 <- paste(genos$bim$V1, genos$bim$V4, genos$bim$V5, genos$bim$V6, "b37", sep = "_")
  genos
}

regularize_ld <- function(ld, ridge = 0.1) {
  ld + ridge * diag(nrow(ld))
}

load_sumstats <- function(opt, chromosomes) {
  model <- opt$model
  stage <- opt$stage
  data_name <- opt$`data-name`

  if (!stage %in% c("01", "12")) stop("--stage must be 01 or 12", call. = FALSE)
  if (!model %in% c("mtag", "gap-trans", "no-prior-lrt", "union+100snps-5e-04-lrt")) {
    stop("Unsupported --model. Use mtag, gap-trans, no-prior-lrt, or union+100snps-5e-04-lrt.", call. = FALSE)
  }

  if (model == "mtag") {
    stop_if_missing_dir(opt$`mtag-dir`, "MTAG directory")
    f <- file.path(opt$`mtag-dir`, sprintf("Sumstat_%s_trait_%s.txt", data_name, ifelse(stage == "01", "1", "2")))
    stop_if_missing(f, "MTAG summary-statistics file")
    df <- as.data.frame(fread(f))
    required <- c("CHR", "BP", "A1", "A2", "mtag_z")
    if (!all(required %in% colnames(df))) stop("MTAG file must contain: ", paste(required, collapse = ", "), call. = FALSE)
    df$id <- paste(df$CHR, df$BP, df$A1, df$A2, "b37", sep = "_")
    df <- df[, c("CHR", "BP", "id", "A1", "A2", "mtag_z")]
  } else if (model == "gap-trans") {
    f <- if (stage == "01") opt$`gap-trans-stage01-file` else opt$`gap-trans-stage12-file`
    stop_if_missing(f, sprintf("gap-trans stage %s file", stage))
    df <- as.data.frame(fread(f))
    required <- c("CHR", "POS", "ID", "ALLELE0", "ALLELE1", "Z_META")
    if (!all(required %in% colnames(df))) stop("gap-trans file must contain: ", paste(required, collapse = ", "), call. = FALSE)
    df <- df[, c("CHR", "POS", "ID", "ALLELE0", "ALLELE1", "Z_META")]
  } else {
    if (model == "no-prior-lrt") stop_if_missing_dir(opt$`mixture-model-dir`, "mixture-model directory")
    if (model == "union+100snps-5e-04-lrt") stop_if_missing_dir(opt$`gap-model-dir`, "GAP model directory")

    df <- data.frame()
    for (chr in chromosomes) {
      if (model == "no-prior-lrt") {
        f <- file.path(opt$`mixture-model-dir`, sprintf("GAP_mixture_%s_chr%s_lrt.txt", data_name, chr))
      } else {
        f <- file.path(opt$`gap-model-dir`, sprintf("GAP_%s_union+100snps-5e-04_chr%s_lrt_new.txt", data_name, chr))
      }
      stop_if_missing(f, sprintf("chromosome %s GWAS file", chr))
      df <- rbind(df, fread(f))
    }

    df <- as.data.frame(df)
    required <- c("snp", "zscore_01", "zscore_12")
    if (!all(required %in% colnames(df))) stop("GAP LRT files must contain: ", paste(required, collapse = ", "), call. = FALSE)

    lbd_01 <- median(df$zscore_01^2, na.rm = TRUE) / qchisq(0.5, df = 1, lower.tail = FALSE)
    lbd_12 <- median(df$zscore_12^2, na.rm = TRUE) / qchisq(0.5, df = 1, lower.tail = FALSE)
    df$zscore_01 <- df$zscore_01 / sqrt(lbd_01)
    df$zscore_12 <- df$zscore_12 / sqrt(lbd_12)

    parts <- strsplit(df$snp, "_", fixed = TRUE)
    df$chr <- as.numeric(vapply(parts, `[`, character(1), 1))
    df$pos <- as.numeric(vapply(parts, `[`, character(1), 2))
    df$ref <- vapply(parts, `[`, character(1), 3)
    df$alt <- vapply(parts, `[`, character(1), 4)

    if (stage == "01") df <- df[, c("chr", "pos", "snp", "ref", "alt", "zscore_01")]
    if (stage == "12") df <- df[, c("chr", "pos", "snp", "ref", "alt", "zscore_12")]
  }

  colnames(df) <- c("chr", "pos", "id", "ref", "alt", "Z")
  as.data.frame(df)
}

run_twas_for_chr <- function(chr, gene_overlap, wgtlist, genos, sumstat_chr, min_r2pred = 0.6, max_impute = 0.5, ld_ridge = 0.1) {
  n_gene <- length(gene_overlap)
  out_tbl <- data.frame(
    ID = character(n_gene), CHR = numeric(n_gene), BEST.GWAS.ID = character(n_gene), BEST.GWAS.Z = numeric(n_gene),
    NSNP = numeric(n_gene), NWGT = numeric(n_gene), TWAS.Z = numeric(n_gene), TWAS.P = numeric(n_gene),
    RANGE = character(n_gene), stringsAsFactors = FALSE
  )

  for (w in seq_along(gene_overlap)) {
    gene <- gene_overlap[w]
    wgt_matrix <- wgtlist[wgtlist$gene == gene, , drop = FALSE]
    pos <- as.numeric(vapply(strsplit(wgt_matrix$snp, "_", fixed = TRUE), `[`, character(1), 2))
    snp_range <- paste0(min(pos, na.rm = TRUE), "-", max(pos, na.rm = TRUE))

    m <- match(wgt_matrix$snp, genos$bim$V2)
    m_keep <- !is.na(m)
    wgt_matrix <- wgt_matrix[m_keep, , drop = FALSE]
    cur_fail <- FALSE
    cur_z <- numeric(0)

    if (nrow(wgt_matrix) == 0) {
      cat("WARNING:", gene, "had no SNPs in the LD reference panel\n")
      cur_fail <- TRUE
    } else {
      cur_genos <- scale(genos$bed[, m[m_keep], drop = FALSE])
      cur_bim <- genos$bim[m[m_keep], , drop = FALSE]
      m_sum <- match(cur_bim$V2, sumstat_chr$id)
      cur_z <- sumstat_chr$Z[m_sum]

      if (sum(wgt_matrix$weight != 0) == 0) {
        cat("WARNING:", gene, "had", length(cur_z), "overlapping SNPs, but none with non-zero expression weights\n")
        cur_fail <- TRUE
      }

      if (length(cur_z) == 0) {
        cat("WARNING:", gene, "had no overlapping SNPs\n")
        cur_fail <- TRUE
      } else if (!cur_fail) {
        cur_ld <- t(cur_genos) %*% cur_genos / (nrow(cur_genos) - 1)
        cur_ld <- regularize_ld(cur_ld, ld_ridge)
        cur_miss <- is.na(cur_z)

        if (sum(cur_miss) != 0) {
          if (sum(!cur_miss) == 0) {
            cat("WARNING:", gene, "had no overlapping GWAS Z-scores, skipping this gene\n")
            cur_fail <- TRUE
          } else if (mean(cur_miss) > max_impute) {
            cat("WARNING:", gene, "had", sum(cur_miss), "/", length(cur_miss), "non-overlapping GWAS Z-scores, skipping this gene\n")
            cur_fail <- TRUE
          } else {
            cur_wgt <- cur_ld[cur_miss, !cur_miss, drop = FALSE] %*% solve(cur_ld[!cur_miss, !cur_miss, drop = FALSE] + ld_ridge * diag(sum(!cur_miss)))
            cur_impz <- cur_wgt %*% cur_z[!cur_miss]
            cur_r2pred <- diag(cur_wgt %*% cur_ld[!cur_miss, !cur_miss, drop = FALSE] %*% t(cur_wgt))
            cur_z[cur_miss] <- cur_impz / sqrt(cur_r2pred)

            all_r2pred <- rep(1, length(cur_z))
            all_r2pred[cur_miss] <- cur_r2pred
            if (sum(is.na(all_r2pred)) != 0) {
              cat("WARNING:", gene, "had missing GWAS Z-scores that could not be imputed, skipping this gene\n")
              cur_fail <- TRUE
            } else if (mean(all_r2pred[wgt_matrix$weight != 0]) < min_r2pred) {
              cat("WARNING:", gene, "had mean GWAS Z-score imputation r2 of", mean(all_r2pred[wgt_matrix$weight != 0]), "at expression-weight SNPs, skipping this gene\n")
              cur_fail <- TRUE
            }
          }
        }

        if (!cur_fail) {
          cur_twasz <- as.numeric(wgt_matrix$weight %*% cur_z)
          cur_twasr2pred <- as.numeric(wgt_matrix$weight %*% cur_ld %*% wgt_matrix$weight)
          if (cur_twasr2pred > 0) {
            cur_twas <- cur_twasz / sqrt(cur_twasr2pred)
          } else {
            cur_fail <- TRUE
            cat("WARNING:", gene, "had zero predictive variance, try a different model\n")
          }
        }
      }
    }

    out_tbl$CHR[w] <- chr
    out_tbl$ID[w] <- gene
    out_tbl$NSNP[w] <- nrow(wgtlist[wgtlist$gene == gene, , drop = FALSE])
    out_tbl$RANGE[w] <- snp_range

    topgwas <- if (length(cur_z) > 0) which.max(cur_z^2) else integer(0)
    if (!cur_fail && length(topgwas) != 0 && !is.na(topgwas)) {
      out_tbl$BEST.GWAS.ID[w] <- wgt_matrix[topgwas, "snp"]
      out_tbl$BEST.GWAS.Z[w] <- cur_z[topgwas]
      out_tbl$NWGT[w] <- nrow(wgt_matrix)
      out_tbl$TWAS.Z[w] <- cur_twas
      out_tbl$TWAS.P[w] <- 2 * pnorm(abs(cur_twas), lower.tail = FALSE)
    } else {
      out_tbl$BEST.GWAS.ID[w] <- NA
      out_tbl$BEST.GWAS.Z[w] <- NA
      out_tbl$NWGT[w] <- ifelse(exists("wgt_matrix"), nrow(wgt_matrix), NA)
      out_tbl$TWAS.Z[w] <- NA
      out_tbl$TWAS.P[w] <- NA
    }
  }

  out_tbl
}

main <- function() {
  start <- Sys.time()
  chromosomes <- parse_chromosomes(opt$chromosomes)
  stop_if_missing(opt$`weights-file`, "weights file")
  if (is.null(opt$`ldref-prefix-template`)) stop("--ldref-prefix-template is required", call. = FALSE)
  if (!dir.exists(opt$`output-dir`)) dir.create(opt$`output-dir`, recursive = TRUE)

  weights_name <- opt$`weights-name`
  if (is.null(weights_name) || weights_name == "") weights_name <- tools::file_path_sans_ext(basename(opt$`weights-file`))

  cat("Loading GWAS summary statistics...\n")
  sumstat <- load_sumstats(opt, chromosomes)

  cat("Loading TWAS weights...\n")
  weights <- as.data.frame(fread(opt$`weights-file`))
  required_weight_cols <- c("gene_id", "snp", "weight")
  if (!all(required_weight_cols %in% colnames(weights))) {
    stop("Weights file must contain: ", paste(required_weight_cols, collapse = ", "), call. = FALSE)
  }
  weights$chr <- as.numeric(vapply(strsplit(weights$snp, "_", fixed = TRUE), `[`, character(1), 1))

  out <- data.frame()
  for (chr in chromosomes) {
    cat("Chromosome", chr, "\n")
    sumstat_chr <- sumstat[sumstat$chr == chr, , drop = FALSE]
    if (nrow(sumstat_chr) == 0) {
      cat("  No GWAS variants for chromosome", chr, "\n")
      next
    }

    ld_prefix <- plink_prefix_for_chr(opt$`ldref-prefix-template`, chr)
    genos <- read_plink_custom(ld_prefix, impute = "avg")
    genos <- harmonize_ld_reference(genos)

    sumstat_chr <- sumstat_chr[match(sumstat_chr$id, genos$bim$V2, nomatch = 0) > 0, , drop = FALSE]
    if (nrow(sumstat_chr) == 0) {
      cat("  No GWAS variants overlap LD reference for chromosome", chr, "\n")
      next
    }

    weight_sub <- weights[weights$chr == chr, , drop = FALSE]
    if (nrow(weight_sub) == 0) {
      cat("  No weights for chromosome", chr, "\n")
      next
    }

    keep <- match(weight_sub$snp, genos$bim$V2, nomatch = 0) > 0
    wgtlist <- weight_sub[keep, c("gene_id", "snp", "weight"), drop = FALSE]
    colnames(wgtlist) <- c("gene", "snp", "weight")
    gene_overlap <- unique(wgtlist$gene)

    if (length(gene_overlap) == 0) {
      cat("  No genes have weights overlapping LD reference for chromosome", chr, "\n")
      next
    }

    res <- run_twas_for_chr(
      chr = chr,
      gene_overlap = gene_overlap,
      wgtlist = wgtlist,
      genos = genos,
      sumstat_chr = sumstat_chr,
      min_r2pred = opt$`min-r2pred`,
      max_impute = opt$`max-impute`,
      ld_ridge = opt$`ld-ridge`
    )

    res$heritable <- 0
    if ("heritable" %in% colnames(weight_sub)) {
      res[res$ID %in% weight_sub$gene_id[weight_sub$heritable == 1], "heritable"] <- 1
    }

    print(summary(res$TWAS.P))
    out <- rbind(out, res)
  }

  out_file <- file.path(opt$`output-dir`, sprintf("%s_%s_%s_%s_TWAS.txt", opt$`data-name`, opt$model, opt$stage, weights_name))
  fwrite(out, out_file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  cat("Saved:", out_file, "\n")
  cat("Elapsed:", as.character(Sys.time() - start), "\n")
}

main()
