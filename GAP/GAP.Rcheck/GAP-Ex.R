pkgname <- "GAP"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('GAP')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("GAP_bayesian_lrt")
### * GAP_bayesian_lrt

flush(stderr()); flush(stdout())

### Name: GAP_bayesian_lrt
### Title: GAP Bayesian Likelihood Ratio Test
### Aliases: GAP_bayesian_lrt

### ** Examples

if (FALSE) {
  prior_params <- GAP_bayesian_prior(input_data, alpha = 0.5, p_threshold = 5e-8, random = 100)
  lrt_results <- GAP_bayesian_lrt(input_data, alpha = 0.5, prior = prior_params)
  head(lrt_results)
}




cleanEx()
nameEx("GAP_bayesian_prior")
### * GAP_bayesian_prior

flush(stderr()); flush(stdout())

### Name: GAP_bayesian_prior
### Title: Estimate GAP Bayesian Prior Parameters
### Aliases: GAP_bayesian_prior

### ** Examples

if (FALSE) {
  set.seed(123)
  n_snps <- 1000
  input_data <- data.frame(
    snp = paste0("rs", 1:n_snps),
    chr = sample(1:22, n_snps, replace = TRUE),
    pos = sample(1:100000000, n_snps),
    beta_01 = rnorm(n_snps, 0, 0.1),
    beta_12 = rnorm(n_snps, 0, 0.1),
    beta_cc = rnorm(n_snps, 0, 0.1),
    se_01 = rep(0.05, n_snps),
    se_12 = rep(0.05, n_snps),
    se_cc = rep(0.05, n_snps)
  )

  prior_params <- GAP_bayesian_prior(input_data, alpha = 0.5, p_threshold = 5e-8, random = 100)
  print(prior_params)
}




cleanEx()
nameEx("find_loci")
### * find_loci

flush(stderr()); flush(stdout())

### Name: find_loci
### Title: Find Independent Loci in GWAS Data
### Aliases: find_loci

### ** Examples

# Example GWAS data
gwas_data <- data.frame(
  chr = c(1, 1, 1, 2, 2),
  pos = c(1000000, 1500000, 3000000, 500000, 2000000),
  pval = c(1e-8, 1e-6, 1e-7, 1e-9, 1e-5)
)

# Find independent loci
loci_result <- find_loci(gwas_data)
print(loci_result)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
