# SLINGtools <img src="man/figures/logo.svg" align="right" height="139" />

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

R functions to facilitate managing, processing and evaluating lipidomics data in R.

## Installation

You can install the development version of SLINGtools like so:

``` r
if (!require("remotes")) install.packages("remotes")
remotes::install_github("SLINGhub/SLINGtools")
```

## Example


``` r
library(SLINGtools)

data_file_path <- system.file("extdata", "Testdata_Lipidomics_MHQuant_Detailed.csv", 
                              package = "SLINGtools")

d <- read_MassHunterCSV(data_file_path)
print(d)
```

## Contributor Code of Conduct

  Please note that the SLINGtools project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.
