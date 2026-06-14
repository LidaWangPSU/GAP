## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(GAP)

gap_input <- gap_example_data()

dim(gap_input)
names(gap_input)
head(gap_input)

## -----------------------------------------------------------------------------
gap_example_path()

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

