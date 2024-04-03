
# kinetics

<!-- badges: start -->
<!-- badges: end -->

The goal of kinetics is to provide handy tools to fit biological data to complex kinetics models. In the development version we provide functions to fit biphasic exponential and monophasic exponential decay models. Currently implementing the Levenberg-Marquardt Nonlinear Least-Squares Algorithm from the minpack.lm package, with optimized initial values input, we are able to obtain robust estimates of the biphasic decay model. Our automatic flow also generates estimates of half-life for overall and seperate phases through root-finding methods as well as throwing fitted plots. 

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

