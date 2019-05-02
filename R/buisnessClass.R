##
#
# Date: 2017-08-04,
# TODO: in estimateObject, make ProbabilityEqual estimation optional!
# TODO: clean the file move the geomm, geomc object to theire correct place!
##
library(mgcv)
library(stringr)

##
# Base object for modelling growth,
# Creates the objects needed for estimation of growth in firms.
#
#' data              - should growth data, possible covariates, previous growth data
#' growthIndex       - index in data (name or col number) indicating which data governs growth
#' growthIndexPrev   - index in data (name or col number) indicating previus time points growth
#' exitIndex         - index in data (name or col number) indicating if the company exited that year
##
growthObjInit <- function(data, growthIndex, growthIndexPrev, exitIndex = NULL){


  res <- list(data = data,
              growth = data[,growthIndex],
              growth_lag  = data[,growthIndexPrev],
              growthIndex = growthIndex,
              growthIndexPrev = growthIndexPrev)
  res$exitIndex = NULL
  if(is.null(exitIndex)==FALSE){
    res$exitIndex = exitIndex
    res$exit      = data[, exitIndex]
  }
  res$estimate <- NULL
  class(res) <- "growthObj"
  return(res)
}

##
# takes a growthObj and estimates the parameters
#
#' growthObj    - the growth object
#' formula      - Formula for modeling the increasing propotion:
#'                $ProbAbove - binomial probabilility for increase
#'                $ProbEqual - binomial probabilility for equal
#'                $ProbExit  - binomial probabilility for exiting
#'                $Above     - formula for number of increased
#'                $Below     - formula for number of deacreased
#' dist         - 2x1 distributrion for decline and growth deafult: geommc, geomm
growthEstimate <- function(growthObj,
                           formula = NULL,
                           dist = list(Above = "geomm", Below = "geommc")
                           ){

  #add if missing Above, Below, ProbEqual ProbAbove
  if(is.null(formula))
  {
    formula = list(Above     = "y ~ 1",
                   Below     = "y ~ 1",
                   ProbEqual = "y ~ 1",
                   ProbAbove = "y ~ 1",
                   ProbExit  = "y ~ 1")
  }
  if(is.null(formula$Above)){
    formula$Above <- as.formula("y ~ 1")
  }
  if(is.null(formula$Below)){
    formula$Below <- as.formula("y ~ 1")
  }
  if(is.null(formula$ProbEqual)){
    formula$ProbEqual <- as.formula("y ~ 1")
  }
  if(is.null(formula$ProbExit)){
    formula$ProbExit <- as.formula("y ~ 1")
  }
  if(is.null(formula$ProbAbove)){
    formula$ProbAbove <- as.formula("y ~ 1")
  }


  growthObj$estimate <- list()
  growthObj$estimate$ProbExit <- NULL
  data       <- growthObj$data
  if(is.null(growthObj$exitIndex) == TRUE){

    growth     <- growthObj$growth
    growth_lag <- growthObj$growth_lag

  }else{
    indexExit  <- data$exit == 1
    dataExit   <- growthObj$data
    dataExit$y <- indexExit
    data       <- growthObj$data[indexExit==F,]
    growth     <- growthObj$growth[indexExit==F]
    growth_lag <- growthObj$growth_lag[indexExit==F]
    growthObj$estimate$ProbExit <- gam(formula = as.formula(formula$ProbExit),
                                       data = dataExit,
                                       family = binomial)
  }

  y <- growth - growth_lag
  indexAbove <- y >  0
  indexProbAbove <- (y > 0) & (growth_lag > 1)
  indexEqual <- y == 0
  indexBelow <- y <  0

  ##
  # estimating the probabilility of growth or steadystate
  ##

  dataEqual    <- data
  dataEqual$y  <- indexEqual
  growthObj$estimate$ProbEqual <- gam(formula = as.formula(formula$ProbEqual),
                                     data = dataEqual,
                                     family = binomial)

  dataAbove    <- data[indexEqual == FALSE & growth_lag > 1,]

  dataAbove$y <- y[indexEqual == FALSE & growth_lag > 1] > 0
  growthObj$estimate$ProbAbove <- gam(formula = as.formula(formula$ProbAbove),
                                      data = dataAbove,
                                      family = binomial)
  ##
  # estimating the distribution for growth or decline
  ##
  data_above   <- data[indexAbove, ]
  data_below   <- data[indexBelow, ]
  data_above$y <-    y[indexAbove]
  data_below$y <-  - y[indexBelow]


  growthObj$estimate$Above <- geomm(formula$Above,  data = data_above)
  growthObj$estimate$Below <- geomcm(formula$Below, data = data_below, growthObj$growth_lag[indexBelow] - 1)
  return(growthObj)
}
##
# takes a growthObj splits the data into to prediction set and estimation set
#
#' growthObj - the growth object
#' predPrec  - prediction data precentage of total object
#' index     - index (name or col number) in growthObj$data indicating where todo the crossvalidation over
##
crossvalGrowth <- function(growthObj, predPrec = 0.2, index = NULL){


  if(is.null(index) == TRUE){
    group = 1:dim(growthObj$data)[1]
  }else{
    group = growthObj$data[, index]
  }
  n_sample = floor(length(group) * predPrec)
  if(n_sample == 0)
  {
    print("warning to few elements to generate prediction set")
    return(-1)
  }
  pred_group = sample(group, n_sample)
  reg_group  = setdiff(group, pred_group)
  if(is.null(index) == TRUE){
    index_pred = pred_group
    index_reg  = reg_group
  }else{
    index_pred = growthObj$data[, index]%in%pred_group
    index_reg  = index_pred == FALSE
  }

  growthObj_out <- growthObj
  growthObj_out$data       = growthObj$data[, index_reg]
  growthObj_out$growth     = growthObj$growth[index_reg]
  growthObj_out$growth_lag = growthObj$growth_lag[index_reg]
  growthObj_out$dataPred   = growthObj$data[, index_pred]
  growthObj_out$growth     = growthObj$growth[, index_pred]
  growthObj_out$growth_lag = growthObj$growth_lag[, index_pred]
  if(is.null(growthObj$exitIndex) ==FALSE)
  {
    growthObj_out$exitPred  = growthObj$exit[index_pred]
    growthObj_out$exit      = growthObj$exit[index_reg]
  }

  return(growthObj_out)
}


##
# probability for buisness class object,
# caculates either the density or the probability for x.
# If supplying newdata it is the average over all the obsevations in the data,
# otherwise it is the average over all the data in the object
##
dgrowth <- function(x, growthObj, newdata = NULL)
{
  if(is.null(growthObj$estimate)){
    print('dgrowth require object to run estimate first\n')
    return(-1)
  }

  if(is.null(newdata ))
    newdata = growthObj$data

  Prob0 <- predict(growthObj$estimate$ProbEqual, newdata = newdata, type = 'response')
  ProbA <- predict(growthObj$estimate$ProbAbove, newdata = newdata, type = 'response')
  pA    <- predict(growthObj$estimate$Above, newdata = newdata)
  pB    <- predict(growthObj$estimate$Below, newdata = newdata)
  if(is.null(growthObj$exitIndex) ==FALSE)
    ProbE <- predict(growthObj$estimate$ProbExit, newdata = newdata, type = 'response')


  growth_lag <- newdata[,growthObj$growthIndexPrev]
  pX <- rep(0, length(x))
  for(i in 1:length(x)){
    x_i = x[i]
    index_equal = x_i == growth_lag
    index_above = x_i > growth_lag
    index_below = (x_i < growth_lag) & (x_i > 0) & (growth_lag > 1)
    pX_i              <- rep(1, length(growth_lag))
    if(is.null(growthObj$exitIndex) ==FALSE){
      if(x_i == 0){
        pX_i[index_zero]          = ProbE[index_zero]
      }else{
        pX_i[index_zero == FALSE] = 1 - ProbE[index_zero == FALSE]
      }
    }
    pX_i[index_equal] <- pX_i[index_equal]  * Prob0[index_equal]
    pX_i[index_above] <- pX_i[index_above]  * (1 - Prob0[index_above])
    pX_i[index_above & growth_lag > 1] <-  pX_i[index_above & growth_lag > 1] *  ProbA[index_above & growth_lag > 1]
    pX_i[index_below] <- pX_i[index_below]  * (1 - Prob0[index_below]) * (1 - ProbA[index_below])
    pX_i[index_below] <- pX_i[index_below]  * dgeomc_( - (x_i - growth_lag[index_below]), pB[index_below], growth_lag[index_below]-1)
    pX_i[index_above] <- pX_i[index_above]  * dgeom_(    (x_i - growth_lag[index_above]), pA[index_above])
    pX[i] <- mean(pX_i)
  }
  return(pX)
}

##
# computes the probabililtiy that an residual of type x is observed,
# when the data is lognormal
#
##
dnormal_diff <- function(x, lmObj, growthPrev, newdata = NULL){

  meanX <- predict(lmObj, newdata = newdata)
  sigma <- summary(lmObj)$sigma
  pX = rep(0, length(x))

  pCorr <- 1- plnorm(0.5,meanlog = meanX, sdlog = sigma)
  for(i in 1:length(x)){

    x_i = x[i]
    pX_i  <- rep(0, length(growthPrev))

    e = growthPrev + x_i
    index0  = e>0
    pX_i[index0] <- plnorm(e[index0]+0.5 ,meanlog = meanX[index0], sdlog = sigma)
    pX_i[index0] <- pX_i[index0] - plnorm(e[index0]-0.5 ,meanlog = meanX[index0], sdlog = sigma)
    pX_i <- pX_i / pCorr
    pX[i] <- mean(pX_i)
  }
  return(pX)
}
##
# computes the ratio for the distribution if the model is lognormal
#
# growth/growthprev
#
# x - the range where the emp/empprev can occuer
##
dnormal_ratio <- function(x, lmObj, growthPrev, newdata = NULL){

  meanX <- predict(lmObj, newdata = newdata)
  meanX <- meanX - log(growthPrev) # Y/growthprev ~ lnorm(mu - logrowthprev, sigma)
  sigma <- summary(lmObj)$sigma
  pX = rep(0, length(x))
  for(i in 1:length(x)){
    pX[i] <- mean(dlnorm(x[i], meanlog =meanX, sdlog = sigma))
  }
  return(pX)
}

##
# difference between year and previous given one does not exit at time t
# TODO: might be incorrect due to not correcting for probability of exit
##
dgrowth_diff <- function(x, growthObj, newdata = NULL)
{
  if(is.null(growthObj$estimate)){
    print('dgrowth require object to run estimate first\n')
    return(-1)
  }

  if(is.null(newdata ))
    newdata = growthObj$data

  Prob0 <- predict(growthObj$estimate$ProbEqual, newdata = newdata, type = 'response')
  ProbA <- predict(growthObj$estimate$ProbAbove, newdata = newdata, type = 'response')
  pA    <- 1/predict(growthObj$estimate$Above, newdata = newdata)
  pB    <- 1/predict(growthObj$estimate$Below, newdata = newdata)
  if(is.null(growthObj$exitIndex) ==FALSE){
    ProbE <- predict(growthObj$estimate$ProbExit, newdata = newdata, type = 'response')
  }


  growth_lag <- newdata[,growthObj$growthIndexPrev]
  pX <- rep(0, length(x))
  for(i in 1:length(x)){

    x_i = x[i]
    pX_i              <- rep(0, length(growth_lag))
    if(x_i < 0){
      index = (growth_lag + x_i ) > 0
      pX_i[index] <- (1 - Prob0[index]) * (1 - ProbA[index])
      index2 <- growth_lag > 2
      pX_i[index2] <- pX_i[index2]  * dgeomc_( -rep(x_i,sum(index2)), pB[index2], growth_lag[index2]-1)
    }else if(x_i > 0){
      index <- growth_lag > 1
      pX_i[index] <- c(1 - Prob0[index]) * c(ProbA[index]) * dgeom_(   x_i, pA[index])
      pX_i[index==FALSE] <- c(1 - Prob0[index==F]) * dgeom_(   x_i, pA[index==F])

    }else{
      pX_i <- Prob0
    }
    pX[i] <- mean(pX_i)
  }
  return(pX)
}
##
#' simulation the evolution N "firm observations", the input is the
#' @param growthObj - a growth object estimated parameters
#' @param X         - (N x K) the covariates used to generate the samples
#'
#' @return Y        - (N x 2) [,1] - exit or not
#'                            [,2] - given no exit size of firm
##
rgrowth <- function(growthObj, X){

  Prob0 <- predict(growthObj$estimate$ProbEqual, newdata = X, type = 'response')
  ProbA <- predict(growthObj$estimate$ProbAbove, newdata = X, type = 'response')
  pA    <- 1/predict(growthObj$estimate$Above, newdata = X)
  pB    <- 1/predict(growthObj$estimate$Below, newdata = X)
  if(is.null(growthObj$exitIndex) ==FALSE)
    ProbE <- predict(growthObj$estimate$ProbExit, newdata = X, type = 'response')


  growth_lag <- X[,growthObj$growthIndexPrev]
  n <- length(growth_lag)

  Y <- matrix(NA, nrow = n, ncol = 2)
  if(is.null(growthObj$exitIndex) ==FALSE){
    U_EXIT <- runif(n)
    Y[ ,1] <- U_EXIT < ProbE
  }else{
    Y[, 1] <- 0
  }
  index_E <- Y[,1] == 0
  nums <- 1:n
  nums_noExit <- nums[index_E]
  n_noExit <- sum(index_E)
  index_0 <- runif(n_noExit) < Prob0[nums_noExit]
  Y[nums_noExit[index_0==T],2] = growth_lag[nums_noExit[index_0==T ]]

  nums_noExit_no0     <- nums_noExit[index_0==0]
  n_noExit_no0  <- length(nums_noExit_no0)
  index_A <- runif(n_noExit_no0) < ProbA[nums_noExit_no0]
  index_A[growth_lag[nums_noExit_no0]==1] = T
  n_A <- sum(index_A)
  if(n_A>0)
    Y[nums_noExit_no0[index_A==T],2] = growth_lag[nums_noExit_no0[index_A==T]] + rgeom(n=n_A,
                                             prob = pA[nums_noExit_no0[index_A==T]])+1
  if(n_noExit_no0-n_A>0)
    Y[nums_noExit_no0[index_A==F],2] = growth_lag[nums_noExit_no0[index_A==F]] - rgeomc(n=n_noExit_no0-n_A,
                                              p = pB[nums_noExit_no0[index_A==F]],
                                              K = growth_lag[nums_noExit_no0[index_A==F]]-1)
  return(Y)
}
