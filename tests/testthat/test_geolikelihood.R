##
# 2017-05-18
# developing test for running
#
##

library(cams)
library(testthat)

context("geometric distribution likelihood test")
##
#
# developing test for constrained parameter space, with simulation for
# to verify likelihood (how?)
#
##
test_that("constrinaed geometric likelihood gradient test, second grad test", {
  #basic setup
  set.seed(1)
  n <- 10
  theta <- rep(5, n)
  p  <- 1 / theta
  K <- 4
  Y <- rgeom(n, p)
  geo_obj_const <- ddlgeo_const(Y, theta, K)
  ##
  # numeric differentiation
  ##
  eps <- 10^-5
  theta_peps <- theta+eps
  d_theta <- (lgeo_const(Y, theta_peps, K) - lgeo_const(Y, theta, K))/eps
  theta_meps <- theta-eps
  dd_theta <- (lgeo_const(Y, theta_peps, K) +
               lgeo_const(Y, theta_meps, K) -
               2 * lgeo_const(Y, theta, K))/eps^2


  expect_equal(d_theta, sum(  attr(geo_obj_const,"gradient")), tolerance=1e-4)
  expect_equal(dd_theta, sum( attr(geo_obj_const,"hessian")), tolerance=1e-4)

  geo_obj_const_cpp <- ddlgeo_const_cpp(Y, theta, rep(K, length(Y)))

  expect_equal(geo_obj_const_cpp$dlogl,  attr(geo_obj_const,"gradient"), tolerance=1e-10)
  expect_equal(geo_obj_const_cpp$ddlogl,  attr(geo_obj_const,"hessian"), tolerance=1e-10)
})

test_that(" geometric likelihood, with covariates, constrained, gradient test, second grad test", {

  set.seed(1)
  n <- 50
  X <- cbind(rep(1, n),
             rnorm(n))
  beta <- c(1,0.1)
  theta <- 1 + exp(X%*%beta)
  K <- 5
  p  <- 1 / theta
  Y <- rgeomc(X, beta, K)
  geo_obj <- dd_lgeom_c(beta, Y, X, K)

  ##
  # numeric diff
  ##
  eps <- 10^-5
  d_beta <- rep(0, 2)
  dd_beta <- matrix(0, nrow = 2, ncol = 2)
  for(i in 1:2){
    beta_eps <- beta
    beta_eps[i] <- beta[i] + eps
    geo_obj_eps <- dd_lgeom_c(beta_eps, Y, X, K)
    d_beta[i] <- ( geo_obj_eps[1]-  geo_obj[1])/eps
    dd_beta[i, ] <-  (attr(geo_obj_eps,"gradient")- attr(geo_obj,"gradient"))/eps
  }
  dd_beta <- (t(dd_beta) + dd_beta)/2
  expect_equal(d_beta, c(attr(geo_obj,"gradient")), tolerance=1e-4)
  expect_equal(dd_beta,attr(geo_obj,"hessian"), tolerance=1e-4)

  #comparing C++ to R
  geo_obj_cpp <- dd_lgeomc_cpp(beta, Y, X, rep(0, length(Y)), rep(K, length(Y)))
  expect_equal(geo_obj_cpp$gradient, c(attr(geo_obj,"gradient")), tolerance=1e-10)
  expect_equal(geo_obj_cpp$hessian, dd_beta,attr(geo_obj,"hessian"), tolerance=1e-10)


})


test_that(" geometric likelihood, with covariates, gradient test, second grad test", {

  set.seed(1)
  n <- 10
  X <- cbind(rep(1, n),
             rnorm(n))
  beta <- c(1,0.1)
  theta <- 1 + exp(X%*%beta)
  p  <- 1 / theta
  Y <- rgeom(n, p)

  geo_obj <- dd_lgeom(beta, Y, X)

  ##
  # numeric diff
  ##
  eps <- 10^-5
  d_beta <- rep(0, 2)
  dd_beta <- matrix(0, nrow = 2, ncol = 2)
  for(i in 1:2){
    beta_eps <- beta
    beta_eps[i] <- beta[i] + eps
    geo_obj_eps <- dd_lgeom(beta_eps, Y, X)
    d_beta[i] <- ( geo_obj_eps[1]-  geo_obj[1])/eps
    dd_beta[i, ] <-  (attr(geo_obj_eps,"gradient")- attr(geo_obj,"gradient"))/eps
  }
  dd_beta <- (t(dd_beta) + dd_beta)/2
  expect_equal(d_beta, c(attr(geo_obj,"gradient")), tolerance=1e-4)

  expect_equal(dd_beta,attr(geo_obj,"hessian"), tolerance=1e-4)

  #comparing C++ to R
  geo_obj_cpp <- dd_lgeom_cpp(beta, Y, X, rep(0, length(Y)))
  expect_equal(geo_obj_cpp$gradient, c(attr(geo_obj,"gradient")), tolerance=1e-10)
  expect_equal(geo_obj_cpp$hessian, dd_beta,attr(geo_obj,"hessian"), tolerance=1e-10)


})

test_that("simple geometric likelihood gradient test, second grad test", {
  ##
  # basic setup
  ##
  set.seed(1)
  n <- 10
  theta <- rep(5, n)
  p  <- 1 / theta
  Y <- rgeom(n, p) + 1

  geo_obj <- dd_lgeo(Y, theta)


  ##
  # numeric differentiation
  ##
  eps <- 10^-5
  theta_peps <- theta+eps
  d_theta <- (lgeo(Y, theta_peps) - lgeo(Y, theta))/eps
  theta_meps <- theta-eps
  dd_theta <- (lgeo(Y, theta_peps) + lgeo(Y, theta_meps) -2 *lgeo(Y, theta))/eps^2


  expect_equal(d_theta, sum(attr(geo_obj,"gradient")), tolerance=1e-4)
  expect_equal(dd_theta, sum(attr(geo_obj,"hessian")), tolerance=1e-4)
})

test_that("simple geometric likelihood gradient test, second grad test (Rcpp)", {
  ##
  # basic setup
  ##
  set.seed(1)
  n <- 10
  theta <- rep(5, n)
  p  <- 1 / theta
  Y <- rgeom(n, p) + 1

  geo_obj <- dd_lgeo_cpp(as.matrix(Y), as.matrix(theta))

  ##
  # numeric differentiation
  ##
  eps <- 10^-5
  theta_peps <- theta+eps
  d_theta <- (lgeo(Y, theta_peps) - lgeo(Y, theta))/eps
  theta_meps <- theta-eps
  dd_theta <- (lgeo(Y, theta_peps) + lgeo(Y, theta_meps) -2 *lgeo(Y, theta))/eps^2
  expect_equal(d_theta, sum(geo_obj[['dlogl']]), tolerance=1e-4)
  expect_equal(dd_theta, sum(geo_obj[['ddlogl']]), tolerance=1e-4)
})
