
# cams
R-Package for analysis of firm growth data

To install:
``` r
# install.packages("remotes")
remotes::install_github("JonasWallin/cams")
```



## Small tutorial
This is a small tutorial on how to model firm growth using the model proposed in [1]. 
* simulate data, this data is generated so that has similar charateristics as the data used in [1]. 
* Estimate parameters from the simulated data.
* Draw inference from the estimate.

[1] N. Ahmed, F. Delmar, J. Wallin Modeling New-Firm Growth and Survival Using Event Magnitude Regression.

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
  * $ln_roa         - log ROA  
  * $assp_hist      - aspiration
  * $assp_hist_over  - $assp_hist * ($assp_hist>0)
  * $assp_hist_below - $assp_hist * ($assp_hist<0)
  * $exit            - Did the company exit the market

We then create a truncated dummy variable of previous numer of employees:
``` r
data <- generate.data(n,size)
data$emp_prev_trunc <- data$emp_prev
data$emp_prev_trunc[data$emp_prev_trunc > trunc_emp_prev] = trunc_emp_prev
data$emp_prev_trunc <- factor(data$emp_prev_trunc)
``` 

### Estimation

Now we build a `cam` object where we set the following options:
* The probability of growth, decline and exit follows a logistic regression model.
* The growth magnitude follows a geometric distribution (with mean defined by formula Above).
* The decline magnitude follows a truncated geometric distribution (with mean defined by formula Below).

We need to spesifiy which are the target variable `"emp"`, which was the previous years target `"emp_prev"` and the firm exit variable `"exit`.
``` r
cam_obj <- growthObjInit(data, "emp", "emp_prev", exitIndex=  "exit")
formula_list <- list(Above =     " y ~ 1  + ln_roa +  age + assp_hist_over + assp_hist_below +  offset(log(emp_prev))",
                  Below =     " y ~ 1  + ln_roa +  age + assp_hist_over + assp_hist_below +  offset(log(emp_prev))",
                  ProbEqual = " y ~ 1  + ln_roa + assp_hist_over + assp_hist_below + emp_prev_trunc ",
                  ProbAbove = " y ~ 1  + ln_roa + assp_hist_over + assp_hist_below ",
                  ProbExit  = " y ~ 1  + ln_roa + assp_hist_over + assp_hist_below + emp_prev_trunc ")
cam_obj <- growthEstimate(cam_obj, formula = formula_list, dist = list(Above = "geomm", Below = "geommc"))
```                   


### Inference
We first can explore how well the model fits the data. This is done by `diagonstic.geomm` that compares the emperical distribution of the data to the fitted model.
``` r
diagonstic.geomm(cam_obj$estimate$Above, n.groups=20)
```   
The results looks as follows:
``` r
****************************************************** 
Probability         estimated (numerical)
P(Y  =    1 )  =       0.2832  ( 0.2713 )  n =  121 
P(Y  =    2 )  =       0.1960  ( 0.2108 )  n =  94 
P(Y  =    3 )  =       0.1376  ( 0.1592 )  n =  71 
P(Y  =    4 )  =       0.0980  ( 0.0830 )  n =  37 
P(Y  =    5 )  =       0.0706  ( 0.0830 )  n =  37 
P(Y  =    6 )  =       0.0515  ( 0.0448 )  n =  20 
P(Y in (  6 ,  8 ]) =  0.0662  ( 0.0448 )  n =  20 
P(Y in (  8 , 11 ]) =  0.0498  ( 0.0583 )  n =  26 
P(Y in ( 11 , 35 ]) =  0.0465  ( 0.0448 )  n =  20 
``` 
All estimated objects are in `cam_obj$estimate`
* `$ProbExit` is a `mgcv::gam` object for probability of leaving
* `$ProbEqual` is a `mgcv::gam` object for probability of staying at the same number of employees
* `$ProbAbove` is a `mgcv::gam` object for probability of increasing number of employees.
* `$Above` is a `cams` object with effect of magnitude of increasing number of employees.
* `$Below` is a `cams` object with effect of magnitude of decreasing number of employees.

So for instance:
``` r
summary(cam_obj$estimate$ProbExit)
```   
gives:
``` r
Family: binomial 
Link function: logit 

Formula:
y ~ 1 + ln_roa + assp_hist_over + assp_hist_below + emp_prev_trunc

Parametric coefficients:
                  Estimate Std. Error z value Pr(>|z|)    
(Intercept)      -2.232507   0.116279 -19.200  < 2e-16 ***
ln_roa           -3.383268   0.375100  -9.020  < 2e-16 ***
assp_hist_over    0.047998   0.008404   5.711 1.12e-08 ***
assp_hist_below   0.026546   0.013627   1.948   0.0514 .  
emp_prev_trunc11 -0.024894   0.152395  -0.163   0.8702    
emp_prev_trunc12  0.147662   0.154145   0.958   0.3381    
emp_prev_trunc13 -0.144305   0.168041  -0.859   0.3905    
emp_prev_trunc14  0.146366   0.164165   0.892   0.3726    
emp_prev_trunc15 -0.094420   0.114590  -0.824   0.4100    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


R-sq.(adj) =  0.013   Deviance explained = 2.33%
UBRE = -0.41723  Scale est. = 1         n = 10000
```   
<!-- badges: start -->
[![R-CMD-check](https://github.com/JonasWallin/cams/workflows/R-CMD-check/badge.svg)](https://github.com/JonasWallin/cams/actions)
<!-- badges: end -->
