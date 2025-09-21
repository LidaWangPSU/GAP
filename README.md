# GAP (Genetic Association with Progression) Package

GAP provides methods to decompose case–control genetic effects into stage-specific components and test for genetic associations with disease progression (e.g., Healthy → ANA+ and ANA+ → SLE).

## Installation

To use this package, ensure you have the required dependencies installed:

```r
install.packages(c("mvtnorm", "MASS","devtools"))
library(mvtnorm)
library(MASS)
library(devtools)
```
Then you could install GAP from the repository here.

```
devtools::install_github("LidaWangPSU/GAP/GAP")
library(GAP)
```

## Functions

### Main Functions

- `GAP_bayesian_prior()` - Estimate GAP Bayesian prior parameters
- `GAP_bayesian_lrt()` - Perform likelihood ratio tests for progression associations


## Usage Example

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

# View results
head(results)
```

## Data Format
input:
- `snp`: SNP identifier
- `chr`: Chromosome number
- `pos`: Genomic position (in base pairs)
- `beta_01`: Effect size for trait 0 -> 1
- `beta_12`: Effect size for trait 1 -> 2
- `beta_cc`: Effect size for case control study
- `se_01`: Standard error for beta_01
- `se_12`: Standard error for beta_12
- `se_cc`: Standard error for beta_cc

alpha: measuring the logarithm of the ratio of sample size of stage 1 and stage 0 in the datasets(e.g. log10(0.05))

p_threshold: P value thershold used for clumping (e.g. 5e-8)

random: number of random insignificant snps used in estimating prior (e.g. 50, 100, etc.)

## Output

The `GAP_bayesian_lrt()` function returns a data frame with:
- `snp`: SNP identifier
- `beta_01`: Updated effect size for trait 0 -> 1
- `beta_12`: Updated effect size for trait 1 -> 2
- `zscore_01`: Z-score for trait 0 -> 1
- `zscore_12`: Z-score for trait 1 -> 2

## Examples & test data

Example input and output:
https://github.com/LidaWangPSU/GAP/tree/main/example_data

You can reproduce the quick start using those files; the example typically completes in <10 minutes on a standard laptop.

## Contact
Lida Wang [lida.wang.96@gmail.com](lida.wang.96@gmail.com)

