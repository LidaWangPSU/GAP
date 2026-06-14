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
devtools::install_github("LidaWangPSU/GAP")
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
- `find_loci()` identifies approximately independent loci by distance-based clumping.

## Input data format

`GAP_bayesian_prior()` and `GAP_bayesian_lrt()` expect an input data frame with:

- `snp`: SNP identifier
- `chr`: Chromosome number
- `pos`: Genomic position in base pairs
- `beta_01`: Effect size for trait 0 -> 1
- `beta_12`: Effect size for trait 1 -> 2
- `beta_cc`: Effect size for the case-control study
- `se_01`: Standard error for `beta_01`
- `se_12`: Standard error for `beta_12`
- `se_cc`: Standard error for `beta_cc`

Arguments:

- `alpha`: log-scale quantity describing the ratio of sample size between stage-specific datasets, for example `log10(0.05)`
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
chr22 example files.

After installing the package, open the vignette with:

```r
vignette("gap-workflow", package = "GAP")
```

Or browse all available vignettes:

```r
browseVignettes("GAP")
```

The vignette demonstrates how to:

- load the packaged example files
- parse chromosome and position information from the example IDs
- build a valid GAP input table
- run `GAP_bayesian_prior()` on a small working subset

## Packaged example data

The package includes two gzipped example files under `inst/extdata`:

- `example_beta_chr22.txt.gz`
- `example_beta_se_chr22.txt.gz`

You can access them from R with:

```r
gap_example_beta_path()
gap_example_se_path()

example_files <- gap_example_data()
names(example_files)
head(example_files$beta)
head(example_files$se)
```

## Example data repository

The original example files are also available in the GitHub repository:

[https://github.com/LidaWangPSU/GAP/tree/main/example_data](https://github.com/LidaWangPSU/GAP/tree/main/example_data)

## Contact

Lida Wang: [lida.wang.96@gmail.com](mailto:lida.wang.96@gmail.com)
