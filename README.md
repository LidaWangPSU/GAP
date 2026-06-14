# GAP (Genetic Association with Progression) Package

`GAP` provides methods to decompose case-control genetic effects into
stage-specific components and test for genetic associations with disease
progression, such as Healthy -> ANA+ and ANA+ -> SLE.

## Installation

Install the required dependencies first:

```r
install.packages(c("mvtnorm", "MASS", "devtools"))
```

Then install `GAP` from GitHub:

```r
devtools::install_github("LidaWangPSU/GAP", build_vignettes = TRUE, force = TRUE)
library(GAP)
```

You can also install from a local source tarball:

```r
install.packages("GAP_0.1.0.tar.gz", repos = NULL, type = "source")
library(GAP)
```

## Main functions

- `GAP_bayesian_prior()` estimates GAP Bayesian prior parameters.
- `GAP_bayesian_lrt()` performs likelihood ratio tests for progression associations.

## Input data format

`GAP_bayesian_prior()` and `GAP_bayesian_lrt()` expect an input data frame with:

- `snp`: SNP identifier
- `chr`: Chromosome number
- `pos`: Genomic position in base pairs
- `beta_01`: Effect size for progression stage 0 -> 1
- `beta_12`: Effect size for progression stage  1 -> 2
- `beta_cc`: Effect size for the case-control study
- `se_01`: Standard error for `beta_01`
- `se_12`: Standard error for `beta_12`
- `se_cc`: Standard error for `beta_cc`

Arguments:

- `alpha`: log-scale quantity describing the ratio of sample size of stage 1 over stage 0, for example `log10(0.05)`
- `p_threshold`: p-value threshold used for clumping, for example `5e-8`
- `random`: number of random non-significant SNPs used when estimating the prior, for example `50` or `100`

## Usage example

```r
library(GAP)

# input <- your_gwas_data

# Step 1: Estimate prior parameters
prior_params <- GAP_bayesian_prior(
  input = input,
  alpha = log10(0.05),
  p_threshold = 5e-8,
  random = 50
)

# Step 2: Perform likelihood ratio tests
results <- GAP_bayesian_lrt(
  input = input,
  alpha = log10(0.05),
  prior = prior_params
)

head(results)
```

## Output

`GAP_bayesian_lrt()` returns a data frame with:

- `snp`: SNP identifier
- `beta_01`: updated effect size for trait 0 -> 1
- `beta_12`: updated effect size for trait 1 -> 2
- `zscore_01`: z-score for trait 0 -> 1
- `zscore_12`: z-score for trait 1 -> 2

## Vignette

This package includes a vignette showing a full workflow using the packaged
example dataset.

After installing the package, open the vignette with:

```r
vignette("gap-workflow", package = "GAP")
```

Or browse all available vignettes:

```r
browseVignettes("GAP")
```

The vignette demonstrates how to:

- load the packaged example input file
- inspect the required input columns
- run GAP functions

## Packaged example data

The package includes a ready-to-use example input file under `inst/extdata`:

- `example.txt`

You can access them from R with:

```r
gap_example_path()

example_input <- gap_example_data()
head(example_input)
```

## Example data repository

The original example input and output files are also available in the GitHub repository:

[https://github.com/LidaWangPSU/GAP/tree/main/example_data](https://github.com/LidaWangPSU/GAP/tree/main/example_data)

## Contact

Lida Wang: [lida.wang.96@gmail.com](mailto:lida.wang.96@gmail.com)
