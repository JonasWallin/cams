
# cams
R-Package for analysis of firm growth data

To install:
``` r
# install.packages("remotes")
remotes::install_github("JonasWallin/cams")
```



## Small tutorial
Here we show how to:
* simulate data, this data is generated so that has similar charateristics as the model in [1]. 
* Estimate parameters from the simulated data.
* Draw inference from the estimate.

### simulation
Here we generate ten thousand observations for medium sized companies employees
``` r
library(cams)
n <- 10^4
size = "medium"
data <- generate.data(n,size)
```
the data consists of
  * $emp            - how many people was employeed at current year
  * $emp_prev       - previous number of employees.
  * $age            - age of the company.
  * $ln_roa         - log return 
  * $assp_hist      - 
  * $assp_hist_over  - $assp_hist * ($assp_hist>0)
  * $assp_hist_below - $assp_hist * ($assp_hist<0)
  * $exit            - Did the company exit the market


[1] 

<!-- badges: start -->
[![R-CMD-check](https://github.com/JonasWallin/cams/workflows/R-CMD-check/badge.svg)](https://github.com/JonasWallin/cams/actions)
<!-- badges: end -->
