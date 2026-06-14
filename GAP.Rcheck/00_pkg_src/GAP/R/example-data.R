#' Packaged example input file for GAP
#'
#' Returns the path to the packaged example dataset distributed with `GAP`.
#' The file is already formatted for direct use with `GAP_bayesian_prior()` and
#' `GAP_bayesian_lrt()`.
#'
#' @return A length-one character vector.
#' @export
gap_example_path <- function() {
  system.file("extdata", "example.txt", package = "GAP")
}

#' Load packaged example data for GAP
#'
#' Reads the bundled example file and returns a data frame that already matches
#' the required input schema for the main `GAP` functions.
#'
#' @return A data frame with columns `snp`, `chr`, `pos`, `beta_01`,
#'   `beta_12`, `beta_cc`, `se_01`, `se_12`, and `se_cc`.
#' @export
gap_example_data <- function() {
  read.delim(gap_example_path(), check.names = FALSE)
}

#' Deprecated helper for the old paired example-beta file
#'
#' This helper is retained for backward compatibility and now points to the
#' single packaged example input file.
#'
#' @return A length-one character vector.
#' @export
gap_example_beta_path <- function() {
  gap_example_path()
}

#' Deprecated helper for the old paired example-se file
#'
#' This helper is retained for backward compatibility and now points to the
#' single packaged example input file.
#'
#' @return A length-one character vector.
#' @export
gap_example_se_path <- function() {
  gap_example_path()
}
