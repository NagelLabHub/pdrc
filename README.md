
# pdrc (population-based analysis of DNA repair capacity)

<!-- badges: start -->
<!-- badges: end -->

This package is a collection of tools for applying high-throughput DNA repair capacity measurements to population based studies. The tools include data processing and standardization, batch effect correction, regression model fitting, and visualization of results.

## How to use
Tutorials on the R package can be find at the [wiki site](https://github.com/NagelLabHub/pdrc/wiki). 

For interactive use, you can access the [web interface](https://tzhai.shinyapps.io/pdrc/). 

## Installation

You can install the development version of pdrc from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("NagelLabHub/pdrc")
```

## Example

This is a basic example on fitting Comet kinetics data with the package:

``` r
library(pdrc)
## import your data from excel with necessary columns "Sample", and several time point columns with preflix "c_" (e.g., "c_0", "c_15", "c_30", "c_60", "c_120")
result_list <- loop_data(data)
result_df <- convert_results_to_dataframe(result_list)
```

Example files with required format can be found in the [example](https://github.com/NagelLabHub/pdrc/tree/master/example) folder.

## To cite package 'pdrc' in publications use:

   (2025). pdrc: A complex toolkit for processing and analyzing high-throughput DNA repair data. R package version x.x.x. URL: https://github.com/NagelLabHub/pdrc
   
   Additionally, you may cite our manuscript (in preparation): Comprehensive Measurement of Inter-Individual Variation in DNA Repair Capacity in Healthy Individuals. 
   
For LaTeX users, use the following BibTeX entry:

@Manual{pdrc2025,
  title = {pdrc: A Comprehensive Toolkit for Processing and Analyzing High-Throughput DNA Repair Data},
  author = {Ting Zhai, Zachary D Nagel},
  year = {2025},
  note = {R package version x.x.x},
  url = {https://github.com/NagelLabHub/pdrc},
}

## Contributors: 

Primary Developer: Ting Zhai (https://github.com/zhaiting)

Affliation: Nagel Lab, Harvard T.H. Chan School of Public Health

