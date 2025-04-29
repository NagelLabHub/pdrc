
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

## Example

Fitting Comet kinetics data with adaptive model selection:

``` r
library(pdrc)
# Import your data from Excel or CSV
# Required columns: "Sample" and time point columns with prefix "c_" (e.g., "c_0", "c_15", "c_30", "c_60", "c_120")

result_list <- loop_comet_data(data, method = "adaptive")
summary_df <- comet_list_to_df(result_list, method = "adaptive")
```

Example files with required format can be found in the [example](https://github.com/NagelLabHub/pdrc/tree/master/example) folder.

## To cite package 'pdrc' in publications use:

   Ting Zhai, Zachary D Nagel (2025). **pdrc**: A complex toolkit for processing and analyzing high-throughput DNA repair data. R package version 1.0.1. URL: https://github.com/NagelLabHub/pdrc
   
   Alternatively, you may cite our manuscript (in preparation): Comprehensive Measurement of Inter-Individual Variation in DNA Repair Capacity in Healthy Individuals. 
   
For LaTeX users, use the following BibTeX entry:

@Manual{pdrc2025,
  title = {pdrc: A Comprehensive Toolkit for Processing and Analyzing High-Throughput DNA Repair Data},
  author = {Ting Zhai and Zachary D Nagel},
  year = {2025},
  note = {R package version 1.0.1},
  url = {https://github.com/NagelLabHub/pdrc},
}

## Contributors: 

Primary Developer: Ting Zhai.

Advisor: Zachary D Nagel.

Date: 04-28-2025.

