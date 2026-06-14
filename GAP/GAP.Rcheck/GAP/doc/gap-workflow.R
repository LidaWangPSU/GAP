## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(GAP)

example_files <- gap_example_data()
beta <- example_files$beta
se <- example_files$se

dim(beta)
names(beta)[1:5]
head(beta[, 1:4])

## -----------------------------------------------------------------------------
gap_example_beta_path()
gap_example_se_path()

## -----------------------------------------------------------------------------
parse_chr22_id <- function(x) {
  chr_pos <- sub("^.*-([0-9]+):([0-9]+)_.*$", "\\1\t\\2", x)
  parts <- strsplit(chr_pos, "\t", fixed = TRUE)
  data.frame(
    chr = as.integer(vapply(parts, `[`, "", 1)),
    pos = as.integer(vapply(parts, `[`, "", 2))
  )
}

loc <- parse_chr22_id(beta$ID)

gap_input <- data.frame(
  snp = beta$ID,
  chr = loc$chr,
  pos = loc$pos,
  beta_01 = beta$Astrocytes,
  beta_12 = beta$Excitatory.neurons,
  beta_cc = beta$bulk,
  se_01 = se$Astrocytes,
  se_12 = se$Excitatory.neurons,
  se_cc = se$bulk,
  check.names = FALSE
)

head(gap_input)

## -----------------------------------------------------------------------------
set.seed(1)
subset_idx <- sample(seq_len(nrow(gap_input)), 250)
prior_example <- GAP_bayesian_prior(
  input = gap_input[subset_idx, ],
  alpha = 0,
  p_threshold = 1e-4,
  random = 100
)

str(prior_example)

## ----eval = FALSE-------------------------------------------------------------
# lrt_results <- GAP_bayesian_lrt(
#   input = gap_input,
#   alpha = 0,
#   prior = prior_example
# )
# 
# head(lrt_results)

