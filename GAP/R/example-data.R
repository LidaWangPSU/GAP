#' Example chr22 effect-size file bundled with GAP
#'
#' Returns the path to the packaged gzipped effect-size matrix for the chr22
#' example data.
#'
#' @return A length-one character vector.
#' @export
gap_example_beta_path <- function() {
  system.file("extdata", "example_beta_chr22.txt.gz", package = "GAP")
}

#' Example chr22 standard-error file bundled with GAP
#'
#' Returns the path to the packaged gzipped standard-error matrix for the chr22
#' example data.
#'
#' @return A length-one character vector.
#' @export
gap_example_se_path <- function() {
  system.file("extdata", "example_beta_se_chr22.txt.gz", package = "GAP")
}

#' Load packaged chr22 example data
#'
#' Reads the bundled chr22 example effect-size and standard-error files and
#' returns them as a named list.
#'
#' @return A list with entries `beta` and `se`.
#' @export
gap_example_data <- function() {
  list(
    beta = utils::read.delim(gap_example_beta_path(), check.names = FALSE),
    se = utils::read.delim(gap_example_se_path(), check.names = FALSE)
  )
}
