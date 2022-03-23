###
# using continous version
#
#
###
library(formula.tools)
library(mgcv)
##
# beta modeling with constraint, the response Y_i \in [0, K_i]
#
# formula - formula object setting up dependent variables
# data    - the data connected to the formula
# K       - (n x 1) constraint maximum number y can decrease
##
betamc_R <- function(formula, data, K){


  init_list <- init_geommc(formula, data, K)
  y = init_list$y/init_list$K
  X = init_list$X
  K = init_list$K
  offset = init_list$offset
  beta0  = init_list$beta0
  names <- dimnames(beta0)[[1]]
  data$old.respose <- data[,formula.tools::lhs.vars(as.formula(formula))]
  data[,formula.tools::lhs.vars(as.formula(formula))] <- y
  model.fit <- mgcv::gam(as.formula(formula), data= data, family = betar)

  res <- list(data = data,
              formula = formula,
              K = K,
              model.fit = model.fit)

  class(res) <- "betamc"
  return(res)
}

gamma_R <- function(formula, data){

  model.fit <- glm(as.formula(formula), data= data, family = Gamma(link="log"))

  res <- list(data = data,
              formula = formula,
              model.fit = model.fit)

  class(res) <- "gamma"
  return(res)
}

summary.gamma <- function(obj){
  summary(obj$model.fit)
}
summary.betacm <- function(obj){
  summary(obj$model.fit)
}
