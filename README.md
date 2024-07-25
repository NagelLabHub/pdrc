
# pdrc (population-based analysis of DNA repair capacity)

<!-- badges: start -->
<!-- badges: end -->

This package is a collection of tools for applying high-throughput DNA repair capacity measurements to population based studies. The tools include data processing and standardization, batch effect correction, regression model fitting, and visualization of results.

Data from two major types of high-throughput DNA repair measurements can be processed and analyzed: 
1. FM-HCR (Fluorescence Multiplex Host Cell Reactivation), a powerful method to quantify multiple DNA repair pathways in living cells. [PNAS 2014](https://doi.org/10.1073/pnas.1401182111) [Nat Protoc 2021](https://doi.org/10.1038/s41596-021-00577-3)
2. CometChip, a high-throughput comet assay allowing parallel quantification of multiple DNA repair activities. [PNAS 2010](https://doi.org/10.1073/pnas.1004056107)

In the development version we provide functions to standardize data and apply batch correction on FM-HCR data. The batch correction is based on the [ComBat](https://doi.org/10.1093/biostatistics/kxj037) algorithm, which is a popular method for removing batch effects in high-throughput data. 

For Comet data specifically, we include a collection of functions to fit biphasic exponential decay models through robust estimation. Available model fitting algorithms include Levenberg-Marquardt Nonlinear Least-Squares Algorithm from the [minpack.lm package](https://CRAN.R-project.org/package=minpack.lm/), and Bayesian Inference from the [brms package](https://github.com/paul-buerkner/brms/). 

The model fitting process is automated with looping through an list of initial best guesses for the parameters. Our automatic flow also generates estimates of half-life for overall and separate phases through root-finding methods. The model fitting results are visualized with ggplot2, and the goodness of fit is evaluated with R-squared and RSE for nls based algorithm and LOOIC (leave-one-out cross-validation information criteria) for Bayesian algorithm. 

## Installation

You can install the development version of pdrc from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("zhaiting/pdrc")
```

## Example

This is a basic example on fitting Comet kinetics data with the package:

``` r
library(kinetics)
## import your data from excel with necessary columns "Sample", and several time point columns with preflix "c_" (e.g., "c_0", "c_15", "c_30", "c_60", "c_120")
result_list <- loop_data(data)
result_df <- convert_results_to_dataframe(result_list)
```

Example files with required format can be found in the [example](https://github.com/zhaiting/pdrc/tree/master/example) folder.

## To cite package 'pdrc' in publications use:

   (2024). pdrc: A complex toolkit for processing and analyzing high-throughput DNA repair data. R package version x.x.x. URL: https://github.com/zhaiting/pdrc

A BibTeX entry for LaTeX users is:

  @Manual{,
    title = {pdrc: A complex toolkit for processing and analyzing high-throughput DNA repair data},
    author = { },
    year = {2024},
    note = {R package version x.x.x},
    url = {https://github.com/zhaiting/pdrc},
  }
