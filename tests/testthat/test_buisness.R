library(cams)
library(testthat)
context("model object for geometric regression")
test_that("test geometric and constrained geometric_optim", {
set.seed(4)
n <- 300
income <- rnorm(n)
beta <- c(1, 0.2)
theta <- 1 + exp(beta[1] + income  * beta[2] )
p  <- 1 / theta
K <- 10
Y <- rgeom(n, p) + 1
Y_c <- rgeomc(cbind(rep(1,n), income),beta,K)

data = data.frame(emp = Y, income = income)
data_c = data.frame(emp = Y_c, income = income)
output <- geomm( emp ~ 1 + income, data)
output_R <- geomm_R( emp ~ 1 + income, data)
outputc <- geomcm( emp ~ 1 + income, data_c, K)
outputc_R <- geomcm_R( emp ~ 1 + income, data_c, K)


expect_equal(outputc_R$beta, outputc$beta, tolerance=1e-3)
expect_equal(output_R$beta,  output$beta, tolerance=1e-3)
expect_equal(c(outputc_R$beta),  c(beta), tolerance=1e-1)
})
