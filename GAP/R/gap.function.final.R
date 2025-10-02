
# Required packages
# @importFrom mvtnorm dmvnorm
# @importFrom MASS ginv
# @importFrom stats pchisq dnorm cor nlminb

#' Find Independent Loci in GWAS Data
#'
#' This function identifies independent loci from GWAS data by clumping SNPs
#' that are within a specified distance threshold (2Mb by default).
#'
#' @param gwas A data frame containing GWAS results with the following columns:
#'   \itemize{
#'     \item chr: Chromosome number
#'     \item pos: Genomic position (in base pairs)
#'     \item pval: P-value for association
#'   }
#'
#' @return A data frame with the same structure as input \code{gwas}, but with
#'   an additional column \code{top} indicating whether each SNP represents
#'   an independent locus (TRUE) or is clumped with a more significant SNP (FALSE).
#'
#' @details The function works by:
#' \enumerate{
#'   \item Sorting SNPs by p-value within each chromosome
#'   \item Marking the most significant SNP in each chromosome as a top locus
#'   \item For subsequent SNPs, checking if they are within 2Mb of any existing top locus
#'   \item If within distance, marking as clumped; otherwise, marking as new top locus
#' }
#'
#' @examples
#' # Example GWAS data
#' gwas_data <- data.frame(
#'   chr = c(1, 1, 1, 2, 2),
#'   pos = c(1000000, 1500000, 3000000, 500000, 2000000),
#'   pval = c(1e-8, 1e-6, 1e-7, 1e-9, 1e-5)
#' )
#'
#' # Find independent loci
#' loci_result <- find_loci(gwas_data)
#' print(loci_result)
#'
#' @export
find_loci <- function(gwas) {
  chrs <- unique(gwas$chr)

  out <- data.frame()
  for (chr in chrs) {
    cat("Check chr ", chr, "\n")
    gwas_sub <- gwas[which(gwas$chr == chr), ]
    gwas_sub <- gwas_sub[order(gwas_sub$pval), ]
    gwas_sub$top <- TRUE

    region <- gwas_sub[1, "pos"]

    if (nrow(gwas_sub) == 1) {
      out <- rbind(out, gwas_sub)
      next
    }
    cat(nrow(gwas_sub), " significant snps in chr ", chr, "\n")
    for (i in 2:nrow(gwas_sub)) {
      tmp <- gwas_sub[i, "pos"]
      dis <- abs(region - tmp)

      if (min(dis) < 2000 * 1000) {
        gwas_sub$top[i] <- F
      } else {
        region <- c(region, gwas_sub[i, "pos"])
      }
    }

    out <- rbind(out, gwas_sub)
  }

  return(out)
}

#' GAP Prior Likelihood Function (Version 2)
#'
#' Calculates the negative log-likelihood for the GAP (Genetic Association with
#' Progression) prior model. This is an internal function used in the optimization
#' process to estimate prior parameters.
#'
#' @param x0 A numeric vector of length 6 containing the parameters to optimize:
#'   \itemize{
#'     \item x0[1]: s1 parameter (standard deviation component 1)
#'     \item x0[2]: s2 parameter (standard deviation component 2)
#'     \item x0[3]: s3 parameter (correlation component)
#'     \item x0[4]: log(p1) - log probability of state 1
#'     \item x0[5]: log(p2) - log probability of state 2
#'     \item x0[6]: log(p3) - log probability of state 3
#'   }
#' @param beta_m A matrix of effect sizes with columns for beta_01, beta_12, beta_cc
#' @param se_m A matrix of standard errors corresponding to beta_m
#' @param alpha A numeric value for the alpha parameter in the GAP model
#' @param rho A numeric value for the correlation parameter
#'
#' @return A numeric value representing the negative log-likelihood
#'
#' @details This function implements the GAP prior model which considers four
#' possible states for genetic associations:
#' \enumerate{
#'   \item State 00: No association with either trait
#'   \item State 01: Association with trait 1 only
#'   \item State 12: Association with trait 2 only
#'   \item State cc: Association with both traits (Progression)
#' }
#'
#' The function calculates the likelihood for each state and combines them
#' using the prior probabilities to obtain the overall likelihood.
#'
#' @seealso \code{\link{GAP_prior_00_integral_v2}}, \code{\link{GAP_prior_01_integral_v2}},
#'   \code{\link{GAP_prior_12_integral_v2}}, \code{\link{GAP_prior_cc_integral_v2}}
#'
#' @export
GAP_prior_v2 <- function(x0, beta_m, se_m, alpha, rho) {
  s1 <- as.numeric(x0[1])
  s2 <- as.numeric(x0[2])
  s3 <- as.numeric(x0[3])

  if (s1 == 0 || is.na(s1)) {
    s1 <- 0.0001
  }

  if (s2 == 0 || is.na(s2)) {
    s2 <- 0.0001
  }

  if (is.na(s3)) {
    s3 <- 0.0001
  }

  m_t <- matrix(c(s1, s3, 0, s2), nrow = 2, ncol = 2)
  m <- t(m_t) %*% m_t
  t1 <- sqrt(m[1, 1])
  t2 <- sqrt(m[2, 2])
  r <- m[1, 2] / (t1 * t2)

  p1 <- exp(x0[4])
  p2 <- exp(x0[5])
  p3 <- exp(x0[6])
  p4 <- 1 - p1 - p2 - p3

  pp <- c(p1, p2, p3, p4)

  integral_00 <- GAP_prior_00_integral_v2(beta = beta_m[, 2:4], sd = se_m[, 2:4], rho, t1, t2, r, alpha)
  integral_01 <- GAP_prior_01_integral_v2(beta = beta_m[, 2:4], sd = se_m[, 2:4], rho, t1, t2, r, alpha)
  integral_12 <- GAP_prior_12_integral_v2(beta = beta_m[, 2:4], sd = se_m[, 2:4], rho, t1, t2, r, alpha)
  integral_cc <- GAP_prior_cc_integral_v2(beta = beta_m[, 2:4], sd = se_m[, 2:4], rho, t1, t2, r, alpha)

  ll <- cbind(integral_00, integral_01, integral_12, integral_cc)
  colSums(ll)
  ll

  summary(ll)
  t <- apply(ll, 1, which.max)
  table(t)

  p_m <- matrix(pp, nrow = 4, ncol = nrow(beta_m))
  p_m <- t(p_m)
  ll_final <- exp(ll) * p_m
  ll_row <- rowSums(ll_final)

  -(sum(log(ll_row)))
}

#' GAP Prior Integral for State 00 (No Association)
#'
#' Calculates the log-likelihood integral for the GAP model state 00,
#' representing no association with either trait.
#'
#' @param beta A matrix of effect sizes with columns for beta_01, beta_12, beta_cc
#' @param sd A matrix of standard errors corresponding to beta
#' @param rho A numeric value for the correlation parameter
#' @param t1 A numeric value for the first prior standard deviation parameter
#' @param t2 A numeric value for the second prior standard deviation parameter
#' @param r A numeric value for the prior correlation parameter
#' @param alpha A numeric value for the alpha parameter (not used in this state)
#'
#' @return A numeric vector of log-likelihood values for each SNP
#'
#' @details This function calculates the likelihood for the null state where
#' there is no association with either trait. It uses a multivariate normal
#' distribution for the first two traits and a univariate normal for the third.
#'
#' @seealso \code{\link{GAP_prior_v2}}
#'
#' @export
GAP_prior_00_integral_v2 <- function(beta, sd, rho, t1, t2, r, alpha) {
  s1 <- mean(sd[, 1])
  s2 <- mean(sd[, 2])
  sigma <- matrix(c(s1^2, s1 * s2 * rho, s1 * s2 * rho, s2^2), nrow = 2, ncol = 2)
  return(dmvnorm(beta[, 1:2], c(0, 0), sigma, log = T) + dnorm(beta[, 3], 0, sd[, 3], log = T))
}

#' GAP Prior Integral for State 01 (Association with Trait 1 Only)
#'
#' Calculates the log-likelihood integral for the GAP model state 01,
#' representing association with trait 1 only (beta_01 != 0, beta_12 = 0).
#'
#' @param beta A matrix of effect sizes with columns for beta_01, beta_12, beta_cc
#' @param sd A matrix of standard errors corresponding to beta
#' @param rho A numeric value for the correlation parameter
#' @param t1 A numeric value for the first prior standard deviation parameter
#' @param t2 A numeric value for the second prior standard deviation parameter
#' @param r A numeric value for the prior correlation parameter
#' @param alpha A numeric value for the alpha parameter (not used in this state)
#'
#' @return A numeric vector of log-likelihood values for each SNP
#'
#' @details This function calculates the likelihood for the state where there
#' is association with trait 1 only. It integrates over the possible values
#' of beta_01 while constraining beta_12 to be zero.
#'
#' @seealso \code{\link{GAP_prior_v2}}
#'
#' @export
GAP_prior_01_integral_v2 <- function(beta, sd, rho, t1, t2, r, alpha) {
  s1 <- mean(sd[, 1])
  s2 <- mean(sd[, 2])
  sigma <- matrix(c(s1^2, s1 * s2 * rho, s1 * s2 * rho, s2^2), nrow = 2, ncol = 2)

  V1 <- mean(sd[, 1])^2
  V2 <- mean(sd[, 2])^2
  Vcc <- mean(sd[, 3])^2
  A <- 1 / ((1 - rho^2) * V2) + 1 / Vcc + 1 / t1^2
  B <- (V2 * beta[, 1] - rho * beta[, 2] * sqrt(V1 * V2)) / (V1 * V2 * (1 - rho^2)) + beta[, 3] / Vcc
  C <- -(V1 * beta[, 2]^2 + V2 * beta[, 1]^2 - 2 * rho * sqrt(V1 * V2) * beta[, 2] * beta[, 1]) / (2 * V1 * V2 * (1 - rho^2)) - 0.5 * beta[, 3]^2 / Vcc
  s <- 1 / (2 * pi * sqrt(det(sigma))) * 1 / (sqrt(2 * pi * Vcc)) * 1 / (sqrt(2 * pi * t1^2))

  return(log(s) + log(sqrt(2 * pi / A)) + B^2 / (2 * A) + C)
}

#' GAP Prior Integral for State 12 (Association with Trait 2 Only)
#'
#' Calculates the log-likelihood integral for the GAP model state 12,
#' representing association with trait 2 only (beta_01 = 0, beta_12 != 0).
#'
#' @param beta A matrix of effect sizes with columns for beta_01, beta_12, beta_cc
#' @param sd A matrix of standard errors corresponding to beta
#' @param rho A numeric value for the correlation parameter
#' @param t1 A numeric value for the first prior standard deviation parameter
#' @param t2 A numeric value for the second prior standard deviation parameter
#' @param r A numeric value for the prior correlation parameter
#' @param alpha A numeric value for the alpha parameter (affects the variance scaling)
#'
#' @return A numeric vector of log-likelihood values for each SNP
#'
#' @details This function calculates the likelihood for the state where there
#' is association with trait 2 only. The alpha parameter affects the variance
#' scaling for the third trait (beta_cc).
#'
#' @seealso \code{\link{GAP_prior_v2}}
#'
#' @export
GAP_prior_12_integral_v2 <- function(beta, sd, rho, t1, t2, r, alpha) {
  s1 <- mean(sd[, 1])
  s2 <- mean(sd[, 2])
  sigma <- matrix(c(s1^2, s1 * s2 * rho, s1 * s2 * rho, s2^2), nrow = 2, ncol = 2)

  V1 <- mean(sd[, 1])^2
  V2 <- mean(sd[, 2])^2
  Vcc <- mean(sd[, 3])^2
  A <- 1 / ((1 - rho^2) * V2) + 1 / Vcc / (1 + exp(alpha))^2 + 1 / t1^2
  B <- (V2 * beta[, 1] - rho * beta[, 2] * sqrt(V1 * V2)) / (V1 * V2 * (1 - rho^2)) + beta[, 3] / Vcc / (1 + exp(alpha))
  C <- -(V1 * beta[, 2]^2 + V2 * beta[, 1]^2 - 2 * rho * sqrt(V1 * V2) * beta[, 2] * beta[, 1]) / (2 * V1 * V2 * (1 - rho^2)) - 0.5 * beta[, 3]^2 / Vcc
  s <- 1 / (2 * pi * sqrt(det(sigma))) * 1 / (sqrt(2 * pi * Vcc)) * 1 / (sqrt(2 * pi * t1^2))

  return(log(s) + log(sqrt(2 * pi / A)) + B^2 / (2 * A) + C)
}

#' GAP Prior Integral for State cc (Pleiotropic Association)
#'
#' Calculates the log-likelihood integral for the GAP model state cc,
#' representing pleiotropic association with both traits (beta_01 != 0, beta_12 != 0).
#'
#' @param beta A matrix of effect sizes with columns for beta_01, beta_12, beta_cc
#' @param sd A matrix of standard errors corresponding to beta
#' @param rho A numeric value for the correlation parameter
#' @param t1 A numeric value for the first prior standard deviation parameter
#' @param t2 A numeric value for the second prior standard deviation parameter
#' @param r A numeric value for the prior correlation parameter
#' @param alpha A numeric value for the alpha parameter (affects the variance scaling)
#'
#' @return A numeric vector of log-likelihood values for each SNP
#'
#' @details This function calculates the likelihood for the pleiotropic state where
#' there is association with both traits. The alpha parameter affects the variance
#' scaling and correlation structure for the third trait (beta_cc).
#'
#' @seealso \code{\link{GAP_prior_v2}}
#'
#' @export
GAP_prior_cc_integral_v2 <- function(beta, sd, rho, t1, t2, r, alpha) {
  omega <- matrix(c(t1^2, t1 * t2 * r, t1 * t2 * r, t2^2), nrow = 2, ncol = 2)

  s1 <- mean(sd[, 1])
  s2 <- mean(sd[, 2])
  sigma <- matrix(c(s1^2, s1 * s2 * rho, s1 * s2 * rho, s2^2), nrow = 2, ncol = 2)
  vcc <- 1 / mean(sd[, 3])^2 * matrix(c(1 / (1 + exp(alpha))^2, 1 / (1 + exp(alpha)), 1 / (1 + exp(alpha)), 1), nrow = 2, ncol = 2)

  A <- ginv(omega) + ginv(sigma) + vcc
  b <- as.matrix(beta[, 1:2]) %*% ginv(sigma) + matrix(c((beta[, 3] / mean(sd[, 3])^2 / (1 + exp(alpha))), (beta[, 3] / mean(sd[, 3])^2)), nrow = nrow(beta), ncol = 2)

  sigma_inv <- ginv(sigma)
  c <- -beta[, 3]^2 / (2 * mean(sd[, 3])^2) - 0.5 * (sigma_inv[1, 1] * beta[, 1]^2 + sigma_inv[2, 2] * beta[, 2]^2 + 2 * sigma_inv[1, 2] * beta[, 1] * beta[, 2])
  s <- 1 / sqrt(2 * pi^4 * det(omega) * det(sigma) * mean(sd[, 3])^2 * det(A))

  return(sapply(seq_len(nrow(b)), function(x) {
    log(s) + 0.5 * t(b[x, ]) %*% ginv(A) %*% b[x, ] + c[x]
  }))
}

#' Estimate GAP Bayesian Prior Parameters
#'
#' Estimates the prior parameters for the GAP (Genetic Association with Progression)
#' model using a subset of significant and random SNPs.
#'
#' @param input A data frame containing GWAS results with the following columns:
#'   \itemize{
#'     \item snp: SNP identifier
#'     \item chr: Chromosome number
#'     \item pos: Genomic position (in base pairs)
#'     \item beta_01: Effect size for trait 1
#'     \item beta_12: Effect size for trait 2
#'     \item beta_cc: Effect size for the combined trait
#'     \item se_01: Standard error for beta_01
#'     \item se_12: Standard error for beta_12
#'     \item se_cc: Standard error for beta_cc
#'   }
#' @param alpha A numeric value for the alpha parameter in the GAP model
#' @param p_theshold A numeric value for the p-value threshold to define significant SNPs
#' @param random An integer specifying the number of random SNPs to include in the analysis
#'
#' @return A list containing:
#'   \itemize{
#'     \item tau: A numeric vector of length 3 containing the estimated prior parameters
#'       (t1, t2, r) for the covariance matrix
#'     \item p: A numeric vector of length 4 containing the estimated prior probabilities
#'       for the four states (00, 01, 12, cc)
#'   }
#'
#' @details This function:
#' \enumerate{
#'   \item Calculates p-values for each trait association
#'   \item Estimates the correlation between traits using non-significant SNPs
#'   \item Identifies independent loci for each trait using clumping
#'   \item Selects a subset of significant and random SNPs for parameter estimation
#'   \item Optimizes the GAP prior parameters using maximum likelihood
#' }
#'
#' @examples
#' # Example usage with simulated data
#' set.seed(123)
#' n_snps <- 1000
#' input_data <- data.frame(
#'   snp = paste0("rs", 1:n_snps),
#'   chr = sample(1:22, n_snps, replace = TRUE),
#'   pos = sample(1:100000000, n_snps),
#'   beta_01 = rnorm(n_snps, 0, 0.1),
#'   beta_12 = rnorm(n_snps, 0, 0.1),
#'   beta_cc = rnorm(n_snps, 0, 0.1),
#'   se_01 = rep(0.05, n_snps),
#'   se_12 = rep(0.05, n_snps),
#'   se_cc = rep(0.05, n_snps)
#' )
#'
#' # Estimate prior parameters
#' prior_params <- GAP_bayesian_prior(input_data, alpha = 0.5, p_theshold = 5e-8, random = 100)
#' print(prior_params)
#'
#' @seealso \code{\link{find_loci}}, \code{\link{GAP_prior_v2}}
#'
#' @export
GAP_bayesian_prior <- function(input, alpha, p_thershold, random) {
  beta <- as.data.frame(input[, c("snp", "beta_01", "beta_12", "beta_cc")])
  se <- as.data.frame(input[, c("snp", "se_01", "se_12", "se_cc")])

  ## Calculate p values
  input$pval_01 <- pchisq((beta$beta_01 / se$se_01)^2, df = 1, lower.tail = F)
  input$pval_12 <- pchisq((beta$beta_12 / se$se_12)^2, df = 1, lower.tail = F)
  input$pval_cc <- pchisq((beta$beta_cc / se$se_cc)^2, df = 1, lower.tail = F)

  ## Calculate correlation from 01 and 12
  cutoff <- 0.05
  snps <- beta$snp[which(input$pval_01 > cutoff & input$pval_12 > cutoff & input$pval_cc > cutoff)]
  rho <- cor(beta[match(snps, beta$snp), 2:4])[2, 1]

  ## Clump snps
  gwas <- input[which(input$pval_cc < p_thershold), c("snp", "chr", "pos", "pval_cc")]
  colnames(gwas) <- c("snp", "chr", "pos", "pval")
  loci_cc <- find_loci(gwas)

  gwas <- input[which(input$pval_01 < p_thershold), c("snp", "chr", "pos", "pval_01")]
  colnames(gwas) <- c("snp", "chr", "pos", "pval")
  loci_01 <- find_loci(gwas)

  gwas <- input[which(input$pval_12 < p_thershold), c("snp", "chr", "pos", "pval_12")]
  colnames(gwas) <- c("snp", "chr", "pos", "pval")
  loci_12 <- find_loci(gwas)

  sig_snps <- unique(c(loci_cc$snp[which(loci_cc$top)],
                       loci_01$snp[which(loci_01$top)],
                       loci_12$snp[which(loci_12$top)]))
  set.seed(123)
  random_snps <- sample(input$snp, size = random, replace = F)
  idx <- match(union(sig_snps, random_snps), beta$snp)

  beta_m <- beta[idx, ]
  se_m <- se[idx, ]

  x0 <- c(mean(se_m[, 2]), mean(se_m[, 3]), 0, -3, -3, -3)

  s <- suppressWarnings(
    nlminb(start = x0, GAP_prior_v2, beta_m = beta_m, se_m = se_m, alpha = alpha, rho = rho,
           lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf), upper = c(Inf, Inf, Inf, 0, 0, 0))
  )

  s$par

  s1 <- s$par[1]
  s2 <- s$par[2]
  s3 <- s$par[3]
  m_t <- matrix(c(s1, s3, 0, s2), nrow = 2, ncol = 2)
  m <- t(m_t) %*% m_t
  t1 <- m[1, 1]^0.5
  t2 <- m[2, 2]^0.5
  r <- m[1, 2] / (t1 * t2)

  out <- list("tau" = c(t1, t2, r), "p" = c(exp(s$par[4:6]), 1 - sum(exp(s$par[4:6]))))
  out
  return(out)
}

#' GAP Bayesian Log-likelihood for State 2 (Trait 1 Only)
#'
#' Calculates the log-likelihood for the GAP model state where there is
#' association with trait 1 only (beta_01 != 0, beta_12 = 0).
#'
#' @param x0 A numeric value representing the effect size for trait 2 (beta_12)
#' @param beta A numeric vector of length 3 containing the observed effect sizes
#'   (beta_01, beta_12, beta_cc)
#' @param sd A numeric vector of length 3 containing the standard errors
#' @param rho A numeric value for the correlation parameter
#' @param alpha A numeric value for the alpha parameter (not used in this function)
#' @param t1 A numeric value for the first prior standard deviation parameter
#' @param t2 A numeric value for the second prior standard deviation parameter
#' @param r A numeric value for the prior correlation parameter (not used in this function)
#'
#' @return A numeric value representing the negative log-likelihood
#'
#' @details This is an internal function used in the optimization process
#' for the GAP model. It calculates the likelihood for the state where
#' beta_01 is free to vary while beta_12 is constrained to be zero.
#'
#' @seealso \code{\link{GAP_bayesian_lrt}}
#'
#' @export
GAP_bayesian_p2_logl <- function(x0, beta, sd, rho, alpha, t1, t2, r) {
  omega <- matrix(c(sd[1]^2, rho * sd[1] * sd[2],
                    rho * sd[1] * sd[2], sd[2]^2), nrow = 2, ncol = 2)
  beta_cc <- x0
  -(dmvnorm(c(0, x0), mean = beta[1:2], omega, log = T) + dnorm(beta_cc, beta[3], sd[3], log = T) + dnorm(x0, 0, t2, log = T))
}

#' GAP Bayesian Log-likelihood for State 3 (Trait 2 Only)
#'
#' Calculates the log-likelihood for the GAP model state where there is
#' association with trait 2 only (beta_01 = 0, beta_12 != 0).
#'
#' @param x0 A numeric value representing the effect size for trait 1 (beta_01)
#' @param beta A numeric vector of length 3 containing the observed effect sizes
#'   (beta_01, beta_12, beta_cc)
#' @param sd A numeric vector of length 3 containing the standard errors
#' @param rho A numeric value for the correlation parameter
#' @param alpha A numeric value for the alpha parameter (affects beta_cc transformation)
#' @param t1 A numeric value for the first prior standard deviation parameter
#' @param t2 A numeric value for the second prior standard deviation parameter (not used)
#' @param r A numeric value for the prior correlation parameter (not used)
#'
#' @return A numeric value representing the negative log-likelihood
#'
#' @details This is an internal function used in the optimization process
#' for the GAP model. It calculates the likelihood for the state where
#' beta_12 is free to vary while beta_01 is constrained to be zero.
#' The alpha parameter affects the transformation of beta_cc.
#'
#' @seealso \code{\link{GAP_bayesian_lrt}}
#'
#' @export
GAP_bayesian_p3_logl <- function(x0, beta, sd, rho, alpha, t1, t2, r) {
  omega <- matrix(c(sd[1]^2, rho * sd[1] * sd[2],
                    rho * sd[1] * sd[2], sd[2]^2), nrow = 2, ncol = 2)
  beta_cc <- log((exp(alpha) + 1) / (exp(-x0) + exp(alpha)))
  -(dmvnorm(c(x0, 0), mean = beta[1:2], omega, log = T) + dnorm(beta_cc, beta[3], sd[3], log = T) + dnorm(x0, 0, t1, log = T))
}

#' GAP Bayesian Log-likelihood for State 4 (Pleiotropic Association)
#'
#' Calculates the log-likelihood for the GAP model state where there is
#' pleiotropic association with both traits (beta_01 != 0, beta_12 != 0).
#'
#' @param x0 A numeric vector of length 2 containing the effect sizes
#'   (beta_01, beta_12)
#' @param beta A numeric vector of length 3 containing the observed effect sizes
#'   (beta_01, beta_12, beta_cc)
#' @param sd A numeric vector of length 3 containing the standard errors
#' @param rho A numeric value for the correlation parameter
#' @param alpha A numeric value for the alpha parameter (affects beta_cc transformation)
#' @param t1 A numeric value for the first prior standard deviation parameter
#' @param t2 A numeric value for the second prior standard deviation parameter
#' @param r A numeric value for the prior correlation parameter
#'
#' @return A numeric value representing the negative log-likelihood
#'
#' @details This is an internal function used in the optimization process
#' for the GAP model. It calculates the likelihood for the pleiotropic state
#' where both beta_01 and beta_12 are free to vary. The alpha parameter
#' affects the transformation of beta_cc.
#'
#' @seealso \code{\link{GAP_bayesian_lrt}}
#'
#' @export
GAP_bayesian_p4_logl <- function(x0, beta, sd, rho, alpha, t1, t2, r) {
  omega <- matrix(c(sd[1]^2, rho * sd[1] * sd[2],
                    rho * sd[1] * sd[2], sd[2]^2), nrow = 2, ncol = 2)
  omega_p <- matrix(c(t1^2, t1 * t2 * r, t1 * t2 * r, t2^2), nrow = 2, ncol = 2)

  beta_cc <- log((exp(alpha) + 1) / (exp(-x0[1]) + exp(alpha))) + x0[2]

  -dmvnorm(x0, mean = beta[1:2], sigma = omega, log = T) - dnorm(beta_cc, beta[3], sd[3], log = T) - dmvnorm(x0, c(0, 0), omega_p, log = T)
}

#' GAP Bayesian Log-likelihood for State 4 with Trait 2 Null
#'
#' Calculates the log-likelihood for the GAP model state 4 with the constraint
#' that beta_12 = 0 (trait 2 null hypothesis).
#'
#' @param x0 A numeric value representing the effect size for trait 1 (beta_01)
#' @param beta A numeric vector of length 3 containing the observed effect sizes
#'   (beta_01, beta_12, beta_cc)
#' @param sd A numeric vector of length 3 containing the standard errors
#' @param rho A numeric value for the correlation parameter
#' @param alpha A numeric value for the alpha parameter (affects beta_cc transformation)
#' @param t1 A numeric value for the first prior standard deviation parameter
#' @param t2 A numeric value for the second prior standard deviation parameter
#' @param r A numeric value for the prior correlation parameter
#'
#' @return A numeric value representing the negative log-likelihood
#'
#' @details This is an internal function used in the likelihood ratio test
#' for the GAP model. It calculates the likelihood for the pleiotropic state
#' with the constraint that beta_12 = 0.
#'
#' @seealso \code{\link{GAP_bayesian_lrt}}
#'
#' @export
GAP_bayesian_p4_logl_1_2_null <- function(x0, beta, sd, rho, alpha, t1, t2, r) {
  omega <- matrix(c(sd[1]^2, rho * sd[1] * sd[2],
                    rho * sd[1] * sd[2], sd[2]^2), nrow = 2, ncol = 2)
  omega_p <- matrix(c(t1^2, t1 * t2 * r, t1 * t2 * r, t2^2), nrow = 2, ncol = 2)

  beta_cc <- log((exp(alpha) + 1) / (exp(-x0) + exp(alpha)))

  -dmvnorm(c(x0, 0), mean = beta[1:2], sigma = omega, log = T) - dnorm(beta_cc, beta[3], sd[3], log = T) - dmvnorm(c(x0, 0), c(0, 0), omega_p, log = T)
}

#' GAP Bayesian Log-likelihood for State 4 with Trait 1 Null
#'
#' Calculates the log-likelihood for the GAP model state 4 with the constraint
#' that beta_01 = 0 (trait 1 null hypothesis).
#'
#' @param x0 A numeric value representing the effect size for trait 2 (beta_12)
#' @param beta A numeric vector of length 3 containing the observed effect sizes
#'   (beta_01, beta_12, beta_cc)
#' @param sd A numeric vector of length 3 containing the standard errors
#' @param rho A numeric value for the correlation parameter
#' @param alpha A numeric value for the alpha parameter (not used in this function)
#' @param t1 A numeric value for the first prior standard deviation parameter
#' @param t2 A numeric value for the second prior standard deviation parameter
#' @param r A numeric value for the prior correlation parameter
#'
#' @return A numeric value representing the negative log-likelihood
#'
#' @details This is an internal function used in the likelihood ratio test
#' for the GAP model. It calculates the likelihood for the pleiotropic state
#' with the constraint that beta_01 = 0.
#'
#' @seealso \code{\link{GAP_bayesian_lrt}}
#'
#' @export
GAP_bayesian_p4_logl_0_1_null <- function(x0, beta, sd, rho, alpha, t1, t2, r) {
  omega <- matrix(c(sd[1]^2, rho * sd[1] * sd[2],
                    rho * sd[1] * sd[2], sd[2]^2), nrow = 2, ncol = 2)
  omega_p <- matrix(c(t1^2, t1 * t2 * r, t1 * t2 * r, t2^2), nrow = 2, ncol = 2)

  beta_cc <- x0
  -dmvnorm(c(0, x0), mean = beta[1:2], sigma = omega, log = T) - dnorm(beta_cc, beta[3], sd[3], log = T) - dmvnorm(c(0, x0), c(0, 0), omega_p, log = T)
}

#' GAP Bayesian Likelihood Ratio Test
#'
#' Performs likelihood ratio tests for the GAP (Genetic Association with Progression)
#' model to test for associations with each trait individually and pleiotropically.
#'
#' @param input A data frame containing GWAS results with the following columns:
#'   \itemize{
#'     \item snp: SNP identifier
#'     \item beta_01: Effect size for trait 1
#'     \item beta_12: Effect size for trait 2
#'     \item beta_cc: Effect size for the combined trait
#'     \item se_01: Standard error for beta_01
#'     \item se_12: Standard error for beta_12
#'     \item se_cc: Standard error for beta_cc
#'   }
#' @param alpha A numeric value for the alpha parameter in the GAP model
#' @param prior A list containing the estimated prior parameters from \code{\link{GAP_bayesian_prior}}:
#'   \itemize{
#'     \item tau: A numeric vector of length 3 containing (t1, t2, r)
#'     \item p: A numeric vector of length 4 containing prior probabilities for the four states
#'   }
#'
#' @return A data frame with the following columns:
#'   \itemize{
#'     \item snp: SNP identifier
#'     \item beta_01: Updated effect size for trait 1
#'     \item beta_12: Updated effect size for trait 2
#'     \item zscore_01: Z-score for trait 1 association
#'     \item zscore_12: Z-score for trait 2 association
#'   }
#'
#' @details This function performs likelihood ratio tests for each SNP to test:
#' \enumerate{
#'   \item Association with trait 1 only (beta_01 != 0, beta_12 = 0)
#'   \item Association with trait 2 only (beta_01 = 0, beta_12 != 0)
#'   \item Pleiotropic association (beta_01 != 0, beta_12 != 0)
#' }
#'
#' The function uses the estimated prior parameters to calculate posterior
#' probabilities and updated effect sizes for each SNP.
#'
#' @examples
#' # First estimate prior parameters
#' prior_params <- GAP_bayesian_prior(input_data, alpha = 0.5, p_theshold = 5e-8, random = 100)
#'
#' # Then perform likelihood ratio tests
#' lrt_results <- GAP_bayesian_lrt(input_data, alpha = 0.5, prior = prior_params)
#' head(lrt_results)
#'
#' @seealso \code{\link{GAP_bayesian_prior}}
#'
#' @export
GAP_bayesian_lrt <- function(input, alpha, prior) {
  beta_m <- as.data.frame(input[, c("snp", "beta_01", "beta_12", "beta_cc")])
  se_m <- as.data.frame(input[, c("snp", "se_01", "se_12", "se_cc")])

  ## Calculate p values
  input$pval_01 <- pchisq((beta_m$beta_01 / se_m$se_01)^2, df = 1, lower.tail = F)
  input$pval_12 <- pchisq((beta_m$beta_12 / se_m$se_12)^2, df = 1, lower.tail = F)
  input$pval_cc <- pchisq((beta_m$beta_cc / se_m$se_cc)^2, df = 1, lower.tail = F)

  ## Calculate correlation from 01 and 12
  cutoff <- 0.05
  snps <- input$snp[which(input$pval_01 > cutoff & input$pval_12 > cutoff & input$pval_cc > cutoff)]
  rho <- cor(beta_m[match(snps, beta_m$snp), 2:4])[2, 1]

  beta_m <- as.data.frame(beta_m)
  se_m <- as.data.frame(se_m)

  cat(paste0(nrow(beta_m), " rows of GWAS in total \n"))

  t1 <- prior$tau[1]
  t2 <- prior$tau[2]
  r <- prior$tau[3]

  p1 <- prior$p[1]
  p2 <- prior$p[2]
  p3 <- prior$p[3]
  p4 <- prior$p[4]

  snp <- beta_m$snp

  beta_new_00 <- data.frame()
  beta_new_01 <- data.frame()
  beta_new_12 <- data.frame()
  beta_new_cc <- data.frame()
  se_new_12_01 <- c()
  se_new_01_12 <- c()
  se_new_cc_01 <- c()
  se_new_cc_12 <- c()

  seps <- seq(5000, nrow(beta_m), by = 5000)

  for (i in seq_len(nrow(beta_m))) {
    if (!is.na(match(i, seps))) {
      cat(paste0("Finished ", i, " SNPs \n"))
    }
    sd <- c()
    sd[1] <- se_m$se_01[i]
    sd[2] <- se_m$se_12[i]
    sigma <- matrix(c(sd[1]^2, rho * sd[1] * sd[2],
                      rho * sd[1] * sd[2], sd[2]^2), nrow = 2, ncol = 2)

    omega_inv <- 1 / mean(se_m$se_cc)^2 * matrix(c(1 / (1 + exp(alpha))^2, 1 / (1 + exp(alpha)), 1 / (1 + exp(alpha)), 1), nrow = 2, ncol = 2)

    x0 <- as.numeric(beta_m[i, 2:3])
    b <- as.numeric(beta_m[i, 2:4])
    sd <- as.numeric(se_m[i, 2:4])

    ## p1
    beta_new_00 <- rbind(beta_new_00, c(0, 0))

    ## p2 b01=0
    ### full model
    s <- nlminb(start = 0, GAP_bayesian_p2_logl, beta = b, sd = sd, rho = rho, alpha = alpha, t1 = t1, t2 = t2, r = r)
    ### reduced model with beta_12=0
    logli <- s$objective
    logli_12 <- GAP_bayesian_p2_logl(x0 = 0, beta = b, sd = sd, rho = rho, alpha = alpha, t1 = t1, t2 = t2, r = r)
    delta_12 <- 2 * (logli_12 - logli)
    if (delta_12 < 0) {
      delta_12 <- 0
    }
    GAP_zscore_1_2 <- delta_12^0.5 * sign(s$par)
    beta_new_01 <- rbind(beta_new_01, c(0, s$par))
    se_new_01_12[i] <- s$par / GAP_zscore_1_2

    ## p3 b12=0
    ### full model
    s <- nlminb(start = x0[1], GAP_bayesian_p3_logl, beta = b, sd = sd, rho = rho, alpha = alpha, t1 = t1, t2 = t2, r = r)
    ### reduced model with beta_12=0
    logli <- s$objective
    logli_01 <- GAP_bayesian_p3_logl(x0 = 0, beta = b, sd = sd, rho = rho, alpha = alpha, t1 = t1, t2 = t2, r = r)
    delta_01 <- 2 * (logli_01 - logli)
    if (delta_01 < 0) {
      delta_01 <- 0
    }
    GAP_zscore_0_1 <- delta_01^0.5 * sign(s$par)
    beta_new_12 <- rbind(beta_new_12, c(s$par, 0))
    se_new_12_01[i] <- s$par / GAP_zscore_0_1

    # p4
    ### full model
    s <- nlminb(start = x0, GAP_bayesian_p4_logl, beta = b, sd = sd, rho = rho, alpha = alpha, t1 = t1, t2 = t2, r = r)
    ### reduced model with beta_12=0
    s_12 <- nlminb(start = x0[1], GAP_bayesian_p4_logl_1_2_null, beta = b, sd = sd, rho = rho, alpha = alpha, t1 = t1, t2 = t2, r = r)
    ### reduced model with beta_01=0
    s_01 <- nlminb(start = x0[2], GAP_bayesian_p4_logl_0_1_null, beta = b, sd = sd, rho = rho, alpha = alpha, t1 = t1, t2 = t2, r = r)

    logli <- s$objective
    logli_01 <- s_01$objective
    logli_12 <- s_12$objective

    delta_01 <- 2 * (logli_01 - logli)
    delta_12 <- 2 * (logli_12 - logli)

    if (delta_01 < 0) {
      delta_01 <- 0
    }

    if (delta_12 < 0) {
      delta_12 <- 0
    }

    GAP_zscore_0_1 <- delta_01^0.5 * sign(s$par[1])
    GAP_zscore_1_2 <- delta_12^0.5 * sign(s$par[2])

    beta_new_cc <- rbind(beta_new_cc, s$par)
    se_new_cc_01[i] <- s$par[1] / GAP_zscore_0_1
    se_new_cc_12[i] <- s$par[2] / GAP_zscore_1_2
  }

  ## Mix prior
  beta_new <- p1 * beta_new_00 + p2 * beta_new_01 + p3 * beta_new_12 + p4 * beta_new_cc
  se_new_all_01 <- sqrt(p3^2 * se_new_12_01^2 + p4^2 * se_new_cc_01^2)
  se_new_all_12 <- sqrt(p2^2 * se_new_01_12^2 + p4^2 * se_new_cc_12^2)

  zscore_01 <- beta_new[, 1] / se_new_all_01
  zscore_12 <- beta_new[, 2] / se_new_all_12

  summary(beta_m[1:100, 2:3] / se_m[1:100, 2:3])
  summary(zscore_01)
  summary(zscore_12)

  df <- cbind(snp, beta_new, zscore_01, zscore_12)
  colnames(df) <- c("snp", "beta_01", "beta_12", "zscore_01", "zscore_12")
  return(df)
}


