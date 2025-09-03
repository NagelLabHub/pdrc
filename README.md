
# pdrc (population-based analysis of DNA repair capacity)

<!-- badges: start -->
<!-- badges: end -->

This package is a collection of tools for applying high-throughput DNA repair capacity measurements to population based studies. 
It provides standardized pipelines for data processing, batch correction, flexible regression model fitting (NLS and Bayesian), adaptive modeling based on biological prior knowledge, and high-quality visualization of results.


## What's New in Version 1.0.1

- **New functions added**: `fit_adaptive_spike_model()` for adaptive spike/decay model selection.
- **Fully modular** fitting options: non-linear least squares (NLS), Bayesian biphasic decay modeling, and adaptive spike/decay model selection.
- **Parallelization-ready**: new `setup_parallel()` function to accelerate Bayesian fitting automatically.
- **Harmonized outputs** across all modeling methods (fit object, parameter estimates, half-lives, annotated plots).
- **Updated plotting**: annotations now include sample name, LOOIC, half-lives, and peak information if applicable.


## How to use

- Tutorials are available at the [wiki site](https://github.com/NagelLabHub/pdrc/wiki).
- For interactive use, access the [web interface](https://tzhai.shinyapps.io/pdrc/).


## Installation

You can install the development version of **pdrc** from [GitHub](https://github.com/NagelLabHub/pdrc) with:

``` r
# install.packages("devtools")
devtools::install_github("NagelLabHub/pdrc")
```

## Quick example

Example files with required format can be found in the [example](https://github.com/NagelLabHub/pdrc/tree/master/example) folder.

Performing standardization and batch-effect correction:

``` r
library(pdrc)
# Import your data from Excel or CSV - see example data
# standardization 
list1 <- RE_to_zscore(dat1, c("var1", “var2”)) #generate a list of the standardized data and scaling metrics
# batch correction using the ComBat approach
dat2 <- batch_correction(dat1, pheno, covariate=c(“age”, “sex”), batch="batch")
```

Fitting Comet kinetics data with adaptive model selection:

``` r
library(pdrc)
# Import your data from Excel or CSV
# Required columns (see example data): "Sample" and time point columns with prefix "c_" (e.g., "c_0", "c_15", "c_30", "c_60", "c_120")

result_list <- loop_comet_data(data, method = "adaptive")
summary_df <- comet_list_to_df(result_list, method = "adaptive")
```

## To cite package 'pdrc' in publications use:

Zhai T, Mazzucato P, Ricciardi C, Christiani DC, Liang L, Samson LD, Chaim IA, Nagel ZD. Comprehensive Measurement of Inter-Individual Variation in DNA Repair Capacity in Healthy Individuals. medRxiv. doi: [10.1101/2025.06.13.25329369](https://doi.org/10.1101/2025.06.13.25329369).
