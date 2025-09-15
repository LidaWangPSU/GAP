# GAP (Genetic Association with Progression) Package

This R package provides functions for GAP (Genetic Association with Progression)..

## Installation

To use this package, ensure you have the required dependencies installed:

```r
install.packages(c("mvtnorm", "MASS"))
```

## Functions

### Main Functions

- `GAP_bayesian_prior()` - Estimate GAP Bayesian prior parameters
- `GAP_bayesian_lrt()` - Perform likelihood ratio tests for progression associations


## Usage Example

```r
# Load the functions

# Prepare your data (must have columns: snp, chr, pos, beta_01, beta_12, beta_cc, se_01, se_12, se_cc)
# input_data <- your_gwas_data

# Step 1: Estimate prior parameters
prior_params <- GAP_bayesian_prior(
  input = input_data, 
  alpha = 0.5, 
  p_theshold = 0.05, 
  random = 50
)

# Step 2: Perform likelihood ratio tests
results <- GAP_bayesian_lrt(
  input = input_data, 
  alpha = 0.5, 
  prior = prior_params
)

# View results
head(results)
```

## Data Format

The input data frame must contain the following columns:
- `snp`: SNP identifier
- `chr`: Chromosome number
- `pos`: Genomic position (in base pairs)
- `beta_01`: Effect size for trait 1
- `beta_12`: Effect size for trait 2
- `beta_cc`: Effect size for the combined trait
- `se_01`: Standard error for beta_01
- `se_12`: Standard error for beta_12
- `se_cc`: Standard error for beta_cc

## Output

The `GAP_bayesian_lrt()` function returns a data frame with:
- `snp`: SNP identifier
- `beta_01`: Updated effect size for trait 1
- `beta_12`: Updated effect size for trait 2
- `zscore_01`: Z-score for trait 1 association
- `zscore_12`: Z-score for trait 2 association


