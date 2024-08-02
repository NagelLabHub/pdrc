
# pdrc (population-based analysis of DNA repair capacity)

<!-- badges: start -->
<!-- badges: end -->

This package is a collection of tools for applying high-throughput DNA repair capacity measurements to population based studies. The tools include data processing and standardization, batch effect correction, regression model fitting, and visualization of results.

## How to use
Tutorials on the R package can be find at the [wiki site](https://github.com/zhaiting/pdrc/wiki). 

For interactive use, you can access the [web interface](https://tzhai.shinyapps.io/pdrc/). 

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
