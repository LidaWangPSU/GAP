#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# prepare_biovu_gwas_phenotypes.R
#
# GitHub/shared version of the BioVU phenotype-preparation code.
#
# This script removes private paths and cluster-specific objects. It prepares
# reusable phenotype tables for GWAS or GAP analyses from:
#   1. GRID/RUID mapping files
#   2. EHR laboratory files
#   3. ICD diagnosis-code files
#   4. Optional genotyped-subject ID list
#
# Main outputs:
#   - cleaned GRID/RUID map
#   - cleaned lab table
#   - ANA-positive / ANA-negative phenotype tables
#   - selected autoantibody phenotype tables
#   - cleaned diagnosis table
#   - diagnosis-derived autoimmune cohorts
#   - optional merged GWAS phenotype table
#
# Example:
# Rscript prepare_biovu_gwas_phenotypes.R \
#   --grid-map1 data/Li_PheWAS_GRID_MAPPING.csv \
#   --grid-map2 data/Li_PheWAS_2_Map.csv \
#   --lab1 data/LI_PheWAS_LABs.csv \
#   --lab2 data/LI_PheWAS_LABs_20170630.csv \
#   --lab3 data/Li_PheWAS_LABs_extra.csv \
#   --icd1 data/LI_PheWAS_ICD9_codes_20170630.csv \
#   --icd2 data/LI_PheWAS_ICD9_codes.csv \
#   --icd3 data/Li_PheWAS_ICDs_1-5.csv \
#   --genotyped-ids data/genotyped_grid_ids.txt \
#   --output-dir results/phenotypes \
#   --prefix BIOVU
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(lubridate)
  library(stringr)
})

option_list <- list(
  make_option("--grid-map1", type = "character", default = NA,
              help = "First GRID/RUID map file. Required."),
  make_option("--grid-map2", type = "character", default = NA,
              help = "Second GRID/RUID map file. Optional."),
  make_option("--lab1", type = "character", default = NA,
              help = "Lab file with columns RUID, LAB_DATE, LAB_SHORT_NAME, LAB_VALUE, LAB_UNITS. Optional."),
  make_option("--lab2", type = "character", default = NA,
              help = "Lab file with LAB_LONG_NAME/LAB_SHORT_NAME/LAB_VALUE columns. Optional."),
  make_option("--lab3", type = "character", default = NA,
              help = "Additional lab file already keyed by GRID. Optional."),
  make_option("--icd1", type = "character", default = NA,
              help = "First ICD diagnosis file. Optional."),
  make_option("--icd2", type = "character", default = NA,
              help = "Second ICD diagnosis file keyed by RUID; will be mapped to GRID. Optional."),
  make_option("--icd3", type = "character", default = NA,
              help = "Additional ICD diagnosis file keyed by GRID. Optional."),
  make_option("--genotyped-ids", type = "character", default = NA,
              help = "Optional one-column file containing GRID IDs available for GWAS."),
  make_option("--output-dir", type = "character", default = "gwas_phenotype_output",
              help = "Output directory. Default: %default"),
  make_option("--prefix", type = "character", default = "BIOVU",
              help = "Prefix for output files. Default: %default"),
  make_option("--date-format", type = "character", default = "mdy",
              help = "Date parser: mdy, ymd, or dmy. Default: %default"),
  make_option("--write-rds", action = "store_true", default = FALSE,
              help = "Also write .rds files in addition to .tsv files."),
  make_option("--verbose", action = "store_true", default = FALSE,
              help = "Print phenotype counts.")
)

opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Utility functions
# ------------------------------------------------------------------------------

is_provided <- function(x) !is.na(x) && nzchar(x)

check_file <- function(path, label) {
  if (!is_provided(path)) stop(label, " is required but was not provided.", call. = FALSE)
  if (!file.exists(path)) stop(label, " does not exist: ", path, call. = FALSE)
  invisible(TRUE)
}

read_optional <- function(path, ...) {
  if (!is_provided(path)) return(NULL)
  if (!file.exists(path)) stop("File does not exist: ", path, call. = FALSE)
  fread(path, ...)
}

parse_date <- function(x, fmt = opt$date_format) {
  x <- as.character(x)
  if (fmt == "mdy") return(lubridate::mdy(x, quiet = TRUE))
  if (fmt == "ymd") return(lubridate::ymd(x, quiet = TRUE))
  if (fmt == "dmy") return(lubridate::dmy(x, quiet = TRUE))
  stop("--date-format must be one of: mdy, ymd, dmy", call. = FALSE)
}

clean_text <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "\r", ".")
  x <- make.names(x)
  x <- tolower(x)
  x <- gsub(";", "", x)
  x
}

clean_value <- function(x) {
  x <- as.character(x)
  x <- tolower(x)
  x <- gsub(";", "", x)
  trimws(x)
}

safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

write_table <- function(x, name) {
  out_tsv <- file.path(opt$output_dir, paste0(opt$prefix, "_", name, ".tsv"))
  fwrite(as.data.table(x), out_tsv, sep = "\t", quote = FALSE, na = "NA")
  if (opt$write_rds) {
    saveRDS(x, file.path(opt$output_dir, paste0(opt$prefix, "_", name, ".rds")))
  }
  invisible(out_tsv)
}

summarise_patient_dates <- function(x, id_col = "GRID", date_col = "LAB_DATE", prefix = "phenotype") {
  if (is.null(x) || nrow(x) == 0) {
    return(data.frame(
      GRID = character(),
      min_date = as.Date(character()),
      max_date = as.Date(character()),
      n_records = integer()
    ))
  }

  x %>%
    mutate(.date = parse_date(.data[[date_col]])) %>%
    filter(!is.na(.data[[id_col]]), !is.na(.date)) %>%
    group_by(.data[[id_col]]) %>%
    summarise(
      "{paste0('min_', prefix, '_date')}" := min(.date, na.rm = TRUE),
      "{paste0('max_', prefix, '_date')}" := max(.date, na.rm = TRUE),
      "{paste0('n_', prefix)}" := dplyr::n(),
      .groups = "drop"
    ) %>%
    rename(GRID = all_of(id_col)) %>%
    as.data.frame()
}

print_count <- function(x, label) {
  if (opt$verbose) message(label, ": ", nrow(x), " rows; ", length(unique(x$GRID)), " unique GRID IDs")
}

# ------------------------------------------------------------------------------
# GRID/RUID map
# ------------------------------------------------------------------------------

check_file(opt$`grid-map1`, "--grid-map1")

grid_to_ruid_map1 <- fread(opt$`grid-map1`, colClasses = "character")
grid_to_ruid_map2 <- read_optional(opt$`grid-map2`, colClasses = "character")

grid_to_ruid_map <- if (is.null(grid_to_ruid_map2)) {
  grid_to_ruid_map1
} else {
  rbindlist(list(grid_to_ruid_map1, grid_to_ruid_map2), fill = TRUE)
}
grid_to_ruid_map <- unique(as.data.frame(grid_to_ruid_map))

if (!all(c("GRID", "RUID") %in% colnames(grid_to_ruid_map))) {
  stop("GRID/RUID mapping files must contain columns named GRID and RUID.", call. = FALSE)
}

write_table(grid_to_ruid_map, "grid_to_ruid_map")

# ------------------------------------------------------------------------------
# Lab processing
# ------------------------------------------------------------------------------

standardise_lab1 <- function(path, grid_map) {
  if (!is_provided(path)) return(NULL)

  lab <- fread(path, sep = ",", quote = "", colClasses = "character")

  # Original BioVU format had no header in some exports.
  if (!all(c("RUID", "LAB_DATE", "LAB_SHORT_NAME", "LAB_VALUE", "LAB_UNITS") %in% colnames(lab))) {
    if (ncol(lab) >= 5) {
      lab <- lab[, 1:5]
      colnames(lab) <- c("RUID", "LAB_DATE", "LAB_SHORT_NAME", "LAB_VALUE", "LAB_UNITS")
    } else {
      stop("--lab1 must have at least five columns.", call. = FALSE)
    }
  }

  lab <- merge(lab, grid_map[, c("RUID", "GRID")], by = "RUID", all.x = FALSE, all.y = FALSE)
  lab <- lab[, c("GRID", "LAB_DATE", "LAB_SHORT_NAME", "LAB_VALUE", "LAB_UNITS")]
  colnames(lab) <- c("GRID", "LAB_DATE", "LAB_NAME", "LAB_VALUE", "LAB_UNIT")
  lab$LAB_NAME <- clean_text(lab$LAB_NAME)
  lab$LAB_VALUE <- clean_value(lab$LAB_VALUE)
  lab$LAB_SOURCE <- "lab1"
  unique(as.data.frame(lab))
}

standardise_lab2 <- function(path) {
  if (!is_provided(path)) return(NULL)

  lab <- fread(path, sep = ",", quote = "", colClasses = "character")
  if (!"GRID" %in% colnames(lab)) {
    stop("--lab2 must contain GRID. If it is keyed by RUID, map it to GRID before running this shared script.", call. = FALSE)
  }

  if (!"LAB_DATE" %in% colnames(lab)) stop("--lab2 must contain LAB_DATE.", call. = FALSE)
  if (!"LAB_LONG_NAME" %in% colnames(lab) && !"LAB_NAME" %in% colnames(lab)) {
    stop("--lab2 must contain LAB_LONG_NAME or LAB_NAME.", call. = FALSE)
  }

  lab_name_col <- if ("LAB_LONG_NAME" %in% colnames(lab)) "LAB_LONG_NAME" else "LAB_NAME"
  lab_unit_col <- if ("LAB_UNITS" %in% colnames(lab)) "LAB_UNITS" else if ("LAB_UNIT" %in% colnames(lab)) "LAB_UNIT" else NA

  out <- data.frame(
    GRID = lab$GRID,
    LAB_DATE = lab$LAB_DATE,
    LAB_NAME = lab[[lab_name_col]],
    LAB_VALUE = lab$LAB_VALUE,
    LAB_UNIT = if (is.na(lab_unit_col)) NA_character_ else lab[[lab_unit_col]],
    stringsAsFactors = FALSE
  )

  out$LAB_NAME <- clean_text(out$LAB_NAME)
  out$LAB_VALUE <- clean_value(out$LAB_VALUE)
  out$LAB_SOURCE <- "lab2"
  unique(out)
}

standardise_lab3 <- function(path) {
  if (!is_provided(path)) return(NULL)

  lab <- fread(path, sep = ",", quote = "", colClasses = "character")
  if (!"GRID" %in% colnames(lab)) stop("--lab3 must contain GRID.", call. = FALSE)
  if (!"LAB_DATE" %in% colnames(lab)) stop("--lab3 must contain LAB_DATE.", call. = FALSE)

  lab_name_col <- intersect(c("LAB_NAME", "LAB_LONG_NAME", "LAB_SHORT_NAME"), colnames(lab))[1]
  if (is.na(lab_name_col)) stop("--lab3 must contain a lab-name column.", call. = FALSE)

  lab_unit_col <- intersect(c("LAB_UNIT", "LAB_UNITS"), colnames(lab))[1]

  out <- data.frame(
    GRID = lab$GRID,
    LAB_DATE = lab$LAB_DATE,
    LAB_NAME = lab[[lab_name_col]],
    LAB_VALUE = lab$LAB_VALUE,
    LAB_UNIT = if (is.na(lab_unit_col)) NA_character_ else lab[[lab_unit_col]],
    stringsAsFactors = FALSE
  )

  out$LAB_NAME <- clean_text(out$LAB_NAME)
  out$LAB_VALUE <- clean_value(out$LAB_VALUE)
  out$LAB_SOURCE <- "lab3"
  unique(out)
}

lab_list <- list(
  standardise_lab1(opt$lab1, grid_to_ruid_map),
  standardise_lab2(opt$lab2),
  standardise_lab3(opt$lab3)
)
lab_list <- lab_list[!vapply(lab_list, is.null, logical(1))]

labs <- if (length(lab_list) > 0) {
  unique(rbindlist(lab_list, fill = TRUE))
} else {
  data.table()
}

if (nrow(labs) > 0) {
  write_table(labs, "clean_labs")
}

# ------------------------------------------------------------------------------
# Lab-derived phenotypes
# ------------------------------------------------------------------------------

get_lab_hits <- function(labs, rules) {
  if (is.null(labs) || nrow(labs) == 0) return(data.frame())

  hit_list <- lapply(rules, function(rule) {
    x <- labs

    if (!is.null(rule$names)) {
      x <- x[x$LAB_NAME %in% rule$names, ]
    }

    if (!is.null(rule$name_regex)) {
      x <- x[stringr::str_detect(x$LAB_NAME, rule$name_regex), ]
    }

    if (!is.null(rule$value_exact)) {
      x <- x[x$LAB_VALUE %in% rule$value_exact, ]
    }

    if (!is.null(rule$value_regex)) {
      x <- x[stringr::str_detect(x$LAB_VALUE, rule$value_regex), ]
    }

    if (!is.null(rule$value_numeric_gt)) {
      x$.value_num <- safe_numeric(x$LAB_VALUE)
      x <- x[!is.na(x$.value_num) & x$.value_num > rule$value_numeric_gt, ]
      x$.value_num <- NULL
    }

    if (!is.null(rule$value_numeric_gte)) {
      x$.value_num <- safe_numeric(x$LAB_VALUE)
      x <- x[!is.na(x$.value_num) & x$.value_num >= rule$value_numeric_gte, ]
      x$.value_num <- NULL
    }

    x
  })

  unique(as.data.frame(rbindlist(hit_list, fill = TRUE)))
}

if (nrow(labs) > 0) {
  ana_pos_rules <- list(
    list(names = c("antinuclear_antibodies", "anti.nuclear_antibodies", "anti.nuclear_antibodies_screen", "ana", "ana_direct"), value_exact = "positive"),
    list(names = c("ana.t", "ana.titer", "ana.titer.1", "ana.titer.2"), value_regex = "1:80|1:160|1:320|1:640|1:1280|1:2560"),
    list(names = c("antinuclear.ab"), value_exact = "positive")
  )

  ana_neg_rules <- list(
    list(names = c("antinuclear_antibodies", "anti.nuclear_antibodies", "anti.nuclear_antibodies_screen", "ana", "ana_direct"), value_exact = "negative"),
    list(names = c("ana.t", "ana.titer", "ana.titer.1", "ana.titer.2"), value_regex = "1:16|1:40"),
    list(names = c("antinuclear.ab"), value_exact = "negative")
  )

  dsdna_rules <- list(
    list(names = c("antdna", "anti.dna.sle.current", "anti.dna.sle.comment"), value_regex = "positive"),
    list(names = c("dnatit", "anti.dna.titer"), value_regex = "1:80|1:160")
  )

  rf_rules <- list(
    list(names = c("rheumatoid_factor", "rf", "rheumatoid.factor"), value_exact = "positive"),
    list(names = c("rf", "rf_titer", "rheumatoid.factor"), value_regex = "1:80|1:160|1:320|1:640|1:1280|1:2560"),
    list(names = c("r.f", "rheumatoid.factor"), value_numeric_gt = 15),
    list(names = c("r.f", "rheumatoid.factor"), value_exact = ">2000")
  )

  ccp_rules <- list(
    list(names = c("ccpg", "cyclic.citrul.pep.igg"), value_exact = "positive"),
    list(names = c("ccpg", "cyclic.citrul.pep.igg"), value_numeric_gt = 20)
  )

  scl70_rules <- list(
    list(names = c("sclab", "scl.70.autoabs.eia"), value_exact = c("positive", "strong")),
    list(names = c("sclab", "scl.70.autoabs.eia"), value_numeric_gt = 40)
  )

  ena_rules <- list(
    list(names = c("ssaigg", "ssab.b", "ssab.a", "ssab", "ssa60g", "ssa52g", "ssa_igg_ab", "ssagrf",
                   "ssbigg", "ssb.wb", "ab.rnp", "rnpab", "rnpena", "rnp.wb", "rnp.qn", "rnp.ab",
                   "u1rnp.snrnp_igg", "rnp_antibody", "u1_rnp.smrnp_igg_abs.iaa",
                   "sjogren.s_anti.ss.b", "sjogren.s_anti.ss.a",
                   "ssa..ro..ena..ab..igg", "ssb..la..ena..ab..igg", "rnp..ena...ab.igg",
                   "rnp..ena...igg", "ssa.52..ro...ena..ab.igg.ref",
                   "ssa.60..ro...ena..ab.igg.ref", "ssa..ro...ena..ab.igg.ref",
                   "u1.snrnp.autoabs.eia"), value_exact = c("positive", "pos", "detected")),
    list(names = c("ssaigg", "ssa..ro..ena..ab..igg", "ssb..la..ena..ab..igg", "rnp..ena...ab.igg"), value_numeric_gte = 1)
  )

  smith_rules <- list(
    list(names = c("smab", "sm..smith..autoabs.eia"), value_exact = c("positive", "pos", "detected")),
    list(names = c("smab", "sm..smith..autoabs.eia"), value_numeric_gte = 1)
  )

  lab_phenotypes <- list(
    ana_pos = summarise_patient_dates(get_lab_hits(labs, ana_pos_rules), prefix = "ana_pos"),
    ana_neg = summarise_patient_dates(get_lab_hits(labs, ana_neg_rules), prefix = "ana_neg"),
    dsdna_pos = summarise_patient_dates(get_lab_hits(labs, dsdna_rules), prefix = "dsdna_pos"),
    rf_pos = summarise_patient_dates(get_lab_hits(labs, rf_rules), prefix = "rf_pos"),
    ccp_pos = summarise_patient_dates(get_lab_hits(labs, ccp_rules), prefix = "ccp_pos"),
    scl70_pos = summarise_patient_dates(get_lab_hits(labs, scl70_rules), prefix = "scl70_pos"),
    ena_pos = summarise_patient_dates(get_lab_hits(labs, ena_rules), prefix = "ena_pos"),
    smith_pos = summarise_patient_dates(get_lab_hits(labs, smith_rules), prefix = "smith_pos")
  )

  for (nm in names(lab_phenotypes)) {
    write_table(lab_phenotypes[[nm]], nm)
    print_count(lab_phenotypes[[nm]], nm)
  }
} else {
  lab_phenotypes <- list()
}

# ------------------------------------------------------------------------------
# ICD diagnosis processing
# ------------------------------------------------------------------------------

standardise_icd1 <- function(path) {
  if (!is_provided(path)) return(NULL)
  x <- fread(path, sep = ",", quote = "", colClasses = "character")
  if (!"GRID" %in% colnames(x)) stop("--icd1 must contain GRID.", call. = FALSE)
  if (!"ICD_DATE" %in% colnames(x)) stop("--icd1 must contain ICD_DATE.", call. = FALSE)
  if (!"ICD_CODE" %in% colnames(x)) stop("--icd1 must contain ICD_CODE.", call. = FALSE)
  if (!"ICD_DESCRIPTION" %in% colnames(x)) x$ICD_DESCRIPTION <- NA_character_
  x[, c("GRID", "ICD_DATE", "ICD_CODE", "ICD_DESCRIPTION")]
}

standardise_icd2 <- function(path, grid_map) {
  if (!is_provided(path)) return(NULL)
  x <- fread(path, sep = ",", quote = "", colClasses = "character")

  if (!all(c("RUID", "ICD_DATE", "ICD_CODE", "ICD_DESCRIPTION") %in% colnames(x))) {
    if (ncol(x) >= 4) {
      x <- x[, 1:4]
      colnames(x) <- c("RUID", "ICD_DATE", "ICD_CODE", "ICD_DESCRIPTION")
    } else {
      stop("--icd2 must have at least four columns.", call. = FALSE)
    }
  }

  x <- merge(x, grid_map[, c("RUID", "GRID")], by = "RUID", all.x = FALSE, all.y = FALSE)
  x[, c("GRID", "ICD_DATE", "ICD_CODE", "ICD_DESCRIPTION")]
}

standardise_icd3 <- function(path) {
  if (!is_provided(path)) return(NULL)
  x <- fread(path, sep = ",", quote = "", colClasses = "character")
  if (!"GRID" %in% colnames(x)) stop("--icd3 must contain GRID.", call. = FALSE)
  if (!"ICD_DATE" %in% colnames(x)) stop("--icd3 must contain ICD_DATE.", call. = FALSE)
  if (!"ICD_CODE" %in% colnames(x)) stop("--icd3 must contain ICD_CODE.", call. = FALSE)
  if (!"ICD_DESCRIPTION" %in% colnames(x)) x$ICD_DESCRIPTION <- NA_character_
  x[, c("GRID", "ICD_DATE", "ICD_CODE", "ICD_DESCRIPTION")]
}

icd_list <- list(
  standardise_icd1(opt$icd1),
  standardise_icd2(opt$icd2, grid_to_ruid_map),
  standardise_icd3(opt$icd3)
)
icd_list <- icd_list[!vapply(icd_list, is.null, logical(1))]

diagnosis <- if (length(icd_list) > 0) {
  x <- unique(rbindlist(icd_list, fill = TRUE))
  x <- x[!stringr::str_detect(x$ICD_CODE, "\\\\|NA|^$"), ]
  x$ICD_CODE <- trimws(as.character(x$ICD_CODE))
  x$ICD_FLAG <- ifelse(is.na(stringr::str_extract(x$ICD_CODE, "[A-Za-z]")), "ICD-9", "ICD-10")
  x$DATE <- parse_date(x$ICD_DATE)
  x
} else {
  data.table()
}

if (nrow(diagnosis) > 0) {
  write_table(diagnosis, "clean_diagnosis")
}

extract_icd_cohort <- function(diagnosis, icd10 = character(), icd9 = character(), name) {
  if (is.null(diagnosis) || nrow(diagnosis) == 0) {
    return(data.frame(GRID = character()))
  }

  get_hits <- function(codes, flag) {
    if (length(codes) == 0) return(diagnosis[0, ])
    pattern <- paste0("^(", paste(stringr::str_replace_all(codes, fixed("."), "\\\\."), collapse = "|"), ")")
    diagnosis %>%
      filter(ICD_FLAG == flag, stringr::str_detect(ICD_CODE, pattern))
  }

  hits <- bind_rows(get_hits(icd10, "ICD-10"), get_hits(icd9, "ICD-9"))

  if (nrow(hits) == 0) {
    return(data.frame(
      GRID = character(),
      !!paste0("date_", name) := as.Date(character()),
      !!paste0("date_", name, "_max") := as.Date(character()),
      !!paste0("n_", name) := integer()
    ))
  }

  hits %>%
    filter(!is.na(GRID), !is.na(DATE)) %>%
    group_by(GRID) %>%
    summarise(
      "{paste0('date_', name)}" := min(DATE, na.rm = TRUE),
      "{paste0('date_', name, '_max')}" := max(DATE, na.rm = TRUE),
      "{paste0('n_', name)}" := dplyr::n(),
      .groups = "drop"
    ) %>%
    as.data.frame()
}

if (nrow(diagnosis) > 0) {
  diagnosis_phenotypes <- list(
    rd = extract_icd_cohort(
      diagnosis,
      icd10 = c("M05", "M06", "M31", "M33", "M34", "M35", "I73.0", "M36.0"),
      icd9 = c("446", "710.1", "710.2", "710.3", "710.4", "714", "725"),
      name = "rd"
    ),
    rd_narrow = extract_icd_cohort(
      diagnosis,
      icd10 = c("M33", "M34"),
      icd9 = c("710.1", "710.3"),
      name = "rd_narrow"
    ),
    sle = extract_icd_cohort(
      diagnosis,
      icd10 = c("M32.1", "M32.8", "M32.9"),
      icd9 = c("710.0"),
      name = "sle"
    ),
    cutaneous_le = extract_icd_cohort(
      diagnosis,
      icd10 = c("L93", "H01.12"),
      icd9 = c("695.4", "373.34"),
      name = "le"
    ),
    ra = extract_icd_cohort(
      diagnosis,
      icd10 = c("M05", "M06.8", "M06.9", "M06.0"),
      icd9 = c("714.0", "714.1", "714.2", "714.81"),
      name = "ra"
    ),
    ms = extract_icd_cohort(
      diagnosis,
      icd10 = c("G35"),
      icd9 = c("340"),
      name = "ms"
    ),
    cis = extract_icd_cohort(
      diagnosis,
      icd10 = c("G36.8", "G36.9", "G37.8", "G37.9"),
      icd9 = c("341.8", "341.9"),
      name = "cis"
    ),
    ssc = extract_icd_cohort(
      diagnosis,
      icd10 = c("M34"),
      icd9 = c("710.1"),
      name = "ssc"
    ),
    raynaud = extract_icd_cohort(
      diagnosis,
      icd10 = c("I73.0"),
      icd9 = c("443.0"),
      name = "raynaud"
    ),
    dermatomyositis = extract_icd_cohort(
      diagnosis,
      icd10 = c("M33.1", "M33.10", "M33.11", "M33.12", "M33.13", "M33.19"),
      icd9 = c("710.3"),
      name = "dm"
    )
  )

  # Progression-style timing variables used in the original exploratory code.
  if (nrow(diagnosis_phenotypes$cutaneous_le) > 0 && nrow(diagnosis_phenotypes$sle) > 0) {
    diagnosis_phenotypes$cutaneous_le <- diagnosis_phenotypes$cutaneous_le %>%
      left_join(diagnosis_phenotypes$sle, by = "GRID") %>%
      mutate(days_le_to_sle = as.numeric(date_sle - date_le))
  }

  if (nrow(diagnosis_phenotypes$cis) > 0 && nrow(diagnosis_phenotypes$ms) > 0) {
    diagnosis_phenotypes$cis <- diagnosis_phenotypes$cis %>%
      left_join(diagnosis_phenotypes$ms, by = "GRID") %>%
      mutate(days_cis_to_ms = as.numeric(date_ms - date_cis))
  }

  if (nrow(diagnosis_phenotypes$raynaud) > 0 && nrow(diagnosis_phenotypes$ssc) > 0) {
    diagnosis_phenotypes$raynaud <- diagnosis_phenotypes$raynaud %>%
      left_join(diagnosis_phenotypes$ssc, by = "GRID") %>%
      mutate(days_raynaud_to_ssc = as.numeric(date_ssc - date_raynaud))
  }

  for (nm in names(diagnosis_phenotypes)) {
    write_table(diagnosis_phenotypes[[nm]], paste0("diagnosis_", nm))
    print_count(diagnosis_phenotypes[[nm]], paste0("diagnosis_", nm))
  }

  patient_follow_up <- diagnosis %>%
    filter(!is.na(GRID), !is.na(DATE)) %>%
    group_by(GRID) %>%
    summarise(
      first_diagnosis_date = min(DATE, na.rm = TRUE),
      last_diagnosis_date = max(DATE, na.rm = TRUE),
      follow_up_years = as.numeric(last_diagnosis_date - first_diagnosis_date) / 365.25,
      n_diagnosis_records = dplyr::n(),
      .groups = "drop"
    ) %>%
    as.data.frame()

  write_table(patient_follow_up, "patient_follow_up")
} else {
  diagnosis_phenotypes <- list()
  patient_follow_up <- data.frame(GRID = character())
}

# ------------------------------------------------------------------------------
# Optional combined GWAS phenotype table
# ------------------------------------------------------------------------------

all_ids <- unique(c(
  grid_to_ruid_map$GRID,
  if (exists("labs") && nrow(labs) > 0) labs$GRID else character(),
  if (exists("diagnosis") && nrow(diagnosis) > 0) diagnosis$GRID else character()
))
gwas_pheno <- data.frame(GRID = all_ids, stringsAsFactors = FALSE)

if (is_provided(opt$`genotyped-ids`)) {
  check_file(opt$`genotyped-ids`, "--genotyped-ids")
  geno_ids <- fread(opt$`genotyped-ids`, header = FALSE, colClasses = "character")[[1]]
  gwas_pheno <- gwas_pheno %>% filter(GRID %in% geno_ids)
  gwas_pheno$genotyped <- 1L
}

add_binary <- function(pheno, table, colname) {
  if (is.null(table) || nrow(table) == 0) {
    pheno[[colname]] <- 0L
  } else {
    pheno[[colname]] <- as.integer(pheno$GRID %in% table$GRID)
  }
  pheno
}

for (nm in names(lab_phenotypes)) {
  gwas_pheno <- add_binary(gwas_pheno, lab_phenotypes[[nm]], nm)
}

for (nm in names(diagnosis_phenotypes)) {
  gwas_pheno <- add_binary(gwas_pheno, diagnosis_phenotypes[[nm]], paste0("dx_", nm))
}

# Common GAP-style stage labels:
#   stage0_control: no ANA positivity and no SLE diagnosis
#   stage1_ana_pos: ANA+ and no SLE diagnosis
#   stage2_sle: SLE diagnosis
if (all(c("ana_pos", "dx_sle") %in% colnames(gwas_pheno))) {
  gwas_pheno$stage0_control <- as.integer(gwas_pheno$ana_pos == 0 & gwas_pheno$dx_sle == 0)
  gwas_pheno$stage1_ana_pos <- as.integer(gwas_pheno$ana_pos == 1 & gwas_pheno$dx_sle == 0)
  gwas_pheno$stage2_sle <- as.integer(gwas_pheno$dx_sle == 1)
  gwas_pheno$gap_stage <- ifelse(gwas_pheno$stage2_sle == 1, 2L,
                          ifelse(gwas_pheno$stage1_ana_pos == 1, 1L, 0L))
}

if (exists("patient_follow_up") && nrow(patient_follow_up) > 0) {
  gwas_pheno <- gwas_pheno %>% left_join(patient_follow_up, by = "GRID")
}

write_table(gwas_pheno, "gwas_phenotype_table")

if (opt$verbose) {
  message("Wrote outputs to: ", normalizePath(opt$output_dir))
  message("Combined GWAS phenotype table: ", nrow(gwas_pheno), " individuals")
}
