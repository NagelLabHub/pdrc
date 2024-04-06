
# kinetics

<!-- badges: start -->
<!-- badges: end -->

This package provides handy tools to process and analyze DNA repair data from high-throughput experiments.

Data from two major types of DNA repair experiments can be processed and analyzed: 
1. FM-HCR (Fluorescence Multiplex Host Cell Reactivation), a powerful method to quantify multiple DNA repair pathways in living cells. [PNAS 2014](https://doi.org/10.1073/pnas.1401182111) [Nat Protoc 2021](https://doi.org/10.1038/s41596-021-00577-3)
2. CometChip, a high-throughput comet assay allowing parallel quantification of multiple DNA repair activities. [PNAS 2010](https://doi.org/10.1073/pnas.1004056107)

In the development version we provide functions to fit biphasic exponential and monophasic exponential decay models. Available model fitting algorithms include Levenberg-Marquardt Nonlinear Least-Squares Algorithm from the [minpack.lm package](hhttps://CRAN.R-project.org/package=minpack.lm), and Bayesian Inference from the [brms package](https://github.com/paul-buerkner/brms/). 

The model fitting process is automated with looping through an list of initial best guesses for the parameters. Our automatic flow also generates estimates of half-life for overall and separate phases through root-finding methods. The model fitting results are visualized with ggplot2, and the goodness of fit is evaluated with R-squared and RSE for nls based algorithm and WAIC for Bayesian algorithm. 

## Installation

You can install the development version of kinetics from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ZzzPIPI/kinetics")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(kinetics)
## import your data from excel with necessary columns "Sample", "c_0", "c_15", "c_30", "c_60", "c_120"
result_list <- loop_data(data)
result_df <- convert_results_to_dataframe(result_list)
```

## To cite package 'kinetics' in publications use:

   (2024). Kinetics: A complex toolkit for processing and analyzing high-throughput DNA repair data. R package version x.x.x. URL: https://github.com/ZzzPIPI/kinetics

A BibTeX entry for LaTeX users is:

  @Manual{,
    title = {Kinetics: A complex toolkit for processing and analyzing high-throughput DNA repair data},
    author = { },
    year = {2024},
    note = {R package version x.x.x},
    url = {https://github.com/ZzzPIPI/kinetics},
  }
