##
# 2017-09-25
# Basic model for parameteric geometric distribution, using the
# mean as parameter with p^{-1}(\beta) = \theta(\beta) = 1 + \exp(X\beta)
##
#basic object
dgeom_ <- function(x, p){
  return((1 - p)^(x-1) * p)
}



lgeo <- function(Y, theta, log.l = TRUE)
{
    p <- 1 / theta
    lik <- (Y - 1) * log( 1 - p) + log(p)
    if(log.l)
      return(sum(lik))
    else
      return( exp(sum(lik)) )
}
##
# log likelihook
# Y - n x 1 data natural number larger then zero
#
dd_lgeo <- function(Y, theta)
{
  p <- 1 / theta
  p_m1 <- 1 / (theta - 1)
  llik <- (Y - 1)  * log(1 - p) + log(p)
  Y_p_p_m1 <- (Y - 1) * p * p_m1
  dlog <- - p + Y_p_p_m1
  ddlog <- - p * dlog - Y_p_p_m1 * p_m1
  result <- sum(llik)

  if(is.na(result))
    result <- -Inf
  attr(result,"gradient")  <- dlog
  attr(result, "hessian") <-  ddlog
  return(result)
}

dd_lgeom <- function(beta, Y, X, offset = NULL)
{
  if(is.null(offset)){
    theta_m1 <- exp(X %*% beta)
  }else{
    theta_m1 <- exp(X %*% beta + offset)
  }
  theta <- 1 + theta_m1
  likelihood_list <- dd_lgeo(Y, theta)
  dlogl <- t(X)%*%( attr(likelihood_list,"gradient") * theta_m1)
  ddlogl = t(X)%*%diag(c( attr(likelihood_list,"gradient") * theta_m1))%*%X
  ddlogl <- ddlogl + t(X)%*%diag(c( attr(likelihood_list,"hessian") * theta_m1*theta_m1))%*%X
  result <- likelihood_list[1]
  attr(result,"gradient") <- dlogl
  attr(result, "hessian") <- ddlogl
  return(result)
}

##
#
# constrained parameteric geometric distribution, using the
# mean as parameter with p^{-1}(\beta) = \theta(\beta) = 1 + \exp(X\beta)
##


##
# likelihood of geometric distriubtion, where Y can't be greater then K
#
# Y     - (n x 1) integer in [1,K]
# theta - (n x 1) parameter
# K     - (int)   parameter, maximum value of Y (i.e Y \leq K)
##
lgeo_const <- function(Y, theta, K, log.l = TRUE)
{
  p <- 1 / theta
  p_m <- 1 - p
  lik <- (Y - 1) * log(p_m) + log(p)
  lik <- lik - log(1 - p_m^K)
  if(log.l)
    return(sum(lik))
  else
    return( exp(sum(lik)) )
}


dgeomc_ <- function(x, p, K){
  C_0 <- 1 - (1 - p)^K
  res = rep(0, length(p))
  index_0 <- C_0 > 0
  res <- (1 - p)^(x-1) * p/C_0
  res[index_0==FALSE] = 0
  return(res)
}

# K     - (int)   parameter, maximum value of Y (i.e Y \leq K)
ddlgeo_const <- function(Y, theta, K)
{
  p <- 1 / theta
  theta_inv_m1 <- 1. / (theta - 1)
  p_m1 <- 1 - p
  # likelihood
  llik <- (Y - 1)  * log(p_m1) + log(p)

  base_const <- 1 - p_m1^K
  llik_c <-  log(base_const)  #constraint correction
  # gradient
  Y_p_t_m1 <- (Y - 1) * p * theta_inv_m1
  dlog <- - p + Y_p_t_m1
  p_m1_powK1 <- p_m1^(K-1)
  dlog_c <- - K * p^2 * p_m1_powK1 / base_const # grad constraint correct

  ddlog   <- - p * dlog - Y_p_t_m1 * theta_inv_m1
  ddlog_c <-  (- 2* p + p^2* (K-1)/p_m1 - dlog_c) *  dlog_c
  result <- sum(llik) - sum(llik_c)
  attr(result,"gradient")  <- dlog - dlog_c
  attr(result, "hessian")  <- ddlog - ddlog_c

  return(result)
}
##
# object for optimization of contrained geometric distribution
# Y \in [1, K]
# K     - (int or n x 1)   parameter, maximum value of Y (i.e Y \leq K)
##
dd_lgeom_c <- function(beta, Y, X, K,  offset = NULL)
{
  if(is.null(offset)){
    theta_m1 <- exp(X %*% beta)
  }else{
    theta_m1 <- exp(X %*% beta + offset)
  }

  theta <- 1 + theta_m1
  likelihood_list <- ddlgeo_const(Y, theta, K)
  dlogl <- t(X)%*%( attr(likelihood_list,"gradient") * theta_m1)
  ddlogl = t(X)%*%diag(c( attr(likelihood_list,"gradient") * theta_m1))%*%X
  ddlogl <- ddlogl + t(X)%*%diag(c( attr(likelihood_list,"hessian") * theta_m1*theta_m1))%*%X
  result <- likelihood_list[1]
  attr(result,"gradient") <- dlogl
  attr(result, "hessian") <- ddlogl
  return(result)
}


##
# simulate geometric constrained
#
# X    - (n x p) covariates
# beta - (p x 1) coeffients
# K    - (1 x 1 or n x 1) constraint for Y_i
##
rgeomc <- function(X=NULL , beta=NULL, K, offset = NULL, p = NULL, n=NULL){


  if(length(K) == 1)
    K = rep(K, n)
  if(is.null(p)){
    n <- dim(X)[1]
    if(is.null(offset)){
      theta <- 1 + exp(X%*%beta)
    }else{
      theta <- 1 + exp(X%*%beta + offset)
    }
    p  <- 1 / theta
  }
  Y <- rgeom(n, p) + 1

  done <- FALSE
  for(i in 1:1000)
  {
    index <- Y > K
    n_index <- sum(index)
    if(n_index > 0){
      Y[index] <- rgeom(sum(index),p [index]) + 1
    }else{
      done <- TRUE
      break
    }
  }
  if(done == FALSE){
    cat('error simulation failed\n')
  }
  return(Y)
}


##
# geometric modeling
# estimate through likelihood the geometric distribution with link function
# p(x) = \theta^{-1}(x), where \theta(x) = 1 + \exp(x)
#
# formula - formula object setting up dependent variables
# data    - the data connected to the formula
# offset
##
geomm <- function(formula, data){
  return(geomm_cpp(formula, data))
}

##
# internal function for cpp estimation
#
##
geomm_cpp <- function(formula, data){

  init_list <- init_geomm(formula, data)
  y = init_list$y
  X = init_list$X
  offset = init_list$offset
  if(is.null(offset))
    offset <- rep(0, length(y))
  beta0  = init_list$beta0
  names <- dimnames(beta0)[[1]]
  if(min(y) < 1)
  {
    cat('error response most be postivie\n')
    return()
  }

  Res <-optim_geomm(beta0,
                    y,
                    X,
                    offset);
  beta = as.matrix(Res[["beta"]])
  hessian <- dd_lgeom_cpp(beta, y, X, offset, TRUE)[['hessian']]

  rownames(beta)    <-  names
  rownames(hessian) <-  names
  colnames(hessian) <-  names

  res <- list(y = y,
              X = X,
              beta = beta,
              offset = offset,
              hessian  = hessian,
              formula = formula)
  class(res) <- "geomm"
  return(res)
}





geomm_R <- function(formula, data, offset = NULL){

  init_list <- init_geomm(formula, data)
  y = init_list$y
  X = init_list$X
  offset = init_list$offset
  beta0  = init_list$beta0
  names <- dimnames(beta0)[[1]]
  if(min(y) < 1)
  {
    cat('error response most be postivie\n')
    return()
  }


  f <- function(beta, Y, X, offset){
    res <- dd_lgeom(beta, Y, X, offset = offset)
    res_out <- c()
    res_out[1] <- -res[1]
    attr(res_out,"gradient")  <- -attr(res,"gradient")
    attr(res, "hessian")  <- -attr(res, "hessian")
    return(res_out)
  }
  optim_out <- nlm(f, beta0, hessian = T, Y = y, X = X, offset = offset)
  if(optim_out$code > 2)
    cat('WARNING optimzation in geomm might have failed\n')

  beta <- as.matrix(optim_out$estimate)

  hessian <- optim_out$hessian
  rownames(beta)     <-  names
  rownames(hessian)  <-  names
  colnames(hessian)  <-  names


  res <- list(y = y,
              X = X,
              beta = beta,
              hessian  = hessian,
              offset = offset,
              formula = formula)
  class(res) <- "geomm"
  return(res)
}

init_geomm <- function(formula, data){

  mf <- model.frame(formula = formula, data = data)
  y  <- as.matrix(model.extract(mf, "response"))

  offset <- model.offset(mf)
  X  <- as.matrix(model.matrix(attr(mf, "terms"), data = mf))
  # inital guess
  if(is.null(offset)){
    beta0 <- solve(t(X)%*%X, t(X)%*%(log(y) ))
  }else{
    beta0 <- solve(t(X)%*%X, t(X)%*%(log(y) - offset))
  }

  return(list(y = y, mf = mf, offset = offset, X = X, beta0 = beta0))
}


##
# geometric modeling with constraint, the response Y_i \in [1, K_i]
# estimate through likelihood the geometric distribution with link function
# p(x) = \theta^{-1}(x), where \theta(x) = 1 + \exp(x)
#
# formula - formula object setting up dependent variables
# data    - the data connected to the formula
# K       - (n x 1) constraint maximum number y can decrease
##
geomcm <- function(formula, data, K){
  return(geomcm_cpp(formula, data, K))
}

geomcm_R <- function(formula, data, K){


  init_list <- init_geommc(formula, data, K)
  y = init_list$y
  X = init_list$X
  K = init_list$K
  offset = init_list$offset
  beta0  = init_list$beta0
  names <- dimnames(beta0)[[1]]
  if(min(y) < 1)
  {
    cat('error response most be postivie\n')
    return()
  }

  f <- function(beta, Y, X, K, offset){

    res <- dd_lgeom_c(beta, Y, X, K, offset)
    res_out <- c()
    res_out[1] <- -res[1]
    attr(res_out,"gradient")  <- -attr(res,"gradient")
    attr(res_out, "hessian")  <- -attr(res, "hessian")
    return(res_out)
  }

  optim_out <- nlm(f, beta0, hessian=TRUE, Y = y,
                   X = X,
                   K = K,
                   offset = offset ,
                   iterlim= 1000)
  if(optim_out$code > 2)
    cat('WARNING optimzation in geomcm might have failed\n')

  beta <- as.matrix(optim_out$estimate)
  rownames(beta)[[1]] <-  dimnames(beta0)[[1]]
  hessian <- optim_out$hessian
  rownames(hessian)[[1]] <-  dimnames(beta0)[[1]]
  colnames(hessian)[[2]] <-  dimnames(beta0)[[1]]
  res <- list(y = y,
              X = X,
              beta = beta,
              hessian  = hessian,
              offset = offset,
              formula = formula,
              K = K)

  class(res) <- "geomcm"
  return(res)
}

##
# internal function for initalization of geommc
#
##
init_geommc <- function(formula, data, K, collectY = TRUE){

  if(length(K) == 1)
    K <- rep(K, length(data[, 1]))
  data_trunc <- data
  data_drop <- droplevels(data)
  mf <- model.frame(formula = formula, data = data_drop)
  if(collectY)
    y  <- as.matrix(model.extract(mf, "response"))
  else
    y <- NULL

  offset <- model.offset(mf)
  X  <- as.matrix(model.matrix(attr(mf, "terms"), data = mf))
  # inital guess
  if(is.null(offset)){
    beta0 <- solve(t(X)%*%X, t(X)%*%(log(y) ))
  }else{
    beta0 <- solve(t(X)%*%X, t(X)%*%(log(y) - offset))
  }

  return(list(y = y, mf = mf, offset = offset, X = X, beta0 = beta0, K = K))
}


##
# internal function for cpp estimation
#
##
geomcm_cpp <- function(formula, data, K, offset = NULL){

  init_list <- init_geommc(formula, data, K)
  y = init_list$y
  X = init_list$X
  K = init_list$K
  offset = init_list$offset
  if(is.null(offset))
    offset <- rep(0, length(y))
  beta0  = init_list$beta0
  names <- dimnames(beta0)[[1]]
  if(min(y) < 1)
  {
    cat('error response most be postivie\n')
    return()
  }
  Res <-optim_geommc(beta0,
                     y,
                     as.matrix(X),
                     K,
                     offset);
  beta = as.matrix(Res[["beta"]])
  hessian <- dd_lgeomc_cpp(beta, y, as.matrix(X), offset, K, TRUE)[['hessian']]

  rownames(beta)    <-  names
  rownames(hessian) <-  names
  colnames(hessian) <-  names

  res <- list(y = y,
              X = X,
              beta = beta,
              offset = offset,
              hessian  = hessian,
              formula = formula)
  class(res) <- "geomcm"
  return(res)
}




##
#
# Generating summary of the variables
##
summary.geomcm <- function(object)
{
  names <- dimnames(object$X)[[2]]
  ans <- object["formula"]
  class(ans) <- "summary.geomcm"
  ans$coefficients <- matrix(NA, 0L,4L)
  dimnames(ans$coefficients) <- list(NULL, c("Estimate",
                                             "Std. Error", "t value", "Pr(>|t|)"))
  est <- object$beta
  se <- sqrt(diag(solve(-object$hessian)))
  tval <- est/se
  rdf <- dim(object$X)[1] - length(object$beta)
  ans$coefficients <- cbind(est, se, tval, 2 * pt(abs(tval),
                                                  rdf, lower.tail = FALSE))
  dimnames(ans$coefficients) <- list(dimnames(object$beta)[[1]],
                                     c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
  ans
}

summary.geomm <- function(object)
{
  names <- dimnames(object$X)[[2]]
  ans <- object["formula"]
  class(ans) <- "summary.geomm"
  ans$coefficients <- matrix(NA, 0L,4L)
  dimnames(ans$coefficients) <- list(NULL, c("Estimate",
                                             "Std. Error", "t value", "Pr(>|t|)"))
  est <- object$beta
  se <- sqrt(diag(solve(-object$hessian)))
  tval <- est/se
  rdf <- dim(object$X)[1] - length(object$beta)
  ans$coefficients <- cbind(est, se, tval, 2 * pt(abs(tval),
                                                  rdf, lower.tail = FALSE))
  dimnames(ans$coefficients) <- list(dimnames(object$beta)[[1]],
                                     c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
  ans
}






print.summary.geomm <-function(object){
  cat('\nFamily: geometric  distribution\n')
  cat('link: mu(x) = 1 + exp(x)\n')
  cat('formula: ', Reduce(paste,deparse(object$formula)),'\n')
  res <- printCoefmat(object$coefficients,P.valueS = T, has.Pvalue = T)
}
print.summary.geomcm <-function(object){
  cat('\nFamily: geometric constrained distribution\n')
  cat('link: mu(x) = 1 + exp(x)\n')
  cat('formula: ', Reduce(paste,deparse(object$formula)),'\n')
  res <- printCoefmat(object$coefficients,P.valueS = T, has.Pvalue = T)
}

##
# function for visualising the effect of a covariates on the residuals of a
# geomm model. Red circles are the mean residual, the blue is the median and the bars are the [10%,80%] interval
#' @param geomm       - the geometeric distribution class
#' @param covariate   - either name, or vector, containing the covariate we want to examine vs the residuals
#' @param n.groups    - how many groups should we split the covariates when visualsing the data,
#'                      note that this does not create n.groups rather the quantiles are split into n.groups,
#'                      then one takes the unique quantiles
#' @param standardize - if the residuals should be divided by the offset
##
covariateAnalysis.geomm <- function(geommObj, covariate, n.groups = NULL, standardize = FALSE){

  if("character"%in%is(covariate)){
    covNames <- colnames(geommObj$X)

    if(covariate%in%covNames == FALSE)
      stop(paste(covariate, ' is not in the covariates \n', sep=''))
    covName <- covariate
    covariate <- geommObj$X[,covariate]
  }else{
    if(length(covariate) != length(geommObj$y))
      stop('length of covariate not equal length of y \n')
    covName <- colnames(covariate)
    if(is.null(covName) == FALSE)
      covName <- "covariate"
  }

  if(is.null(n.groups)){
    cov.q <- sort( unique(covariate))
  }else{
    cov.q <- unique(quantile(covariate, seq(0, 1, length=n.groups)))
  }
  #colnames(geomm$X)
  predict.y <- predict(geommObj)
  res <- geommObj$y - predict.y
  if(standardize)
    res <- res / exp(geommObj$offset)

  n.quan <- length(cov.q)
  visual_data <- matrix(0, nrow= n.quan, ncol = 5)
  cov_index <- covariate <= cov.q[1]
  res_temp <- res[cov_index]
  res_quan_temp <- quantile(res_temp, c(0.1, 0.5, 0.9))
  visual_data[1, 1] <- median(covariate[cov_index])
  visual_data[1, 2] <- mean(res_temp)
  visual_data[1, 3] <- res_quan_temp[1]
  visual_data[1, 4] <- res_quan_temp[2]
  visual_data[1, 5] <- res_quan_temp[3]
  for(i in 2:n.quan){
    cov_index <- (covariate > cov.q[i-1]) & (covariate <= cov.q[i])
    res_temp <- res[cov_index]
    res_quan_temp <- quantile(res_temp, c(0.1, 0.5, 0.9))
    visual_data[i, 1] <- median(covariate[cov_index])
    visual_data[i, 2] <- mean(res_temp)
    visual_data[i, 3] <- res_quan_temp[1]
    visual_data[i, 4] <- res_quan_temp[2]
    visual_data[i, 5] <- res_quan_temp[3]
  }

  colnames(visual_data)<-c("median_cov",
                           "mean_res",
                           "q1_res",
                           "q2_res",
                           "q3_res")
  visual_data <- as.data.frame(visual_data)
  g_plot <- ggplot2::ggplot(visual_data, ggplot2::aes(x= median_cov, y=mean_res))  +
            ggplot2::geom_errorbar(ggplot2::aes(ymin = q1_res, ymax = q3_res, color="CI80"),colour="black", width=.1) +
            ggplot2::geom_point(color="red") +
            ggplot2::geom_point(ggplot2::aes(x = median_cov, y = q2_res), shape=5, color="blue") +
            ggplot2::ylab("residual") +
            ggplot2::xlab(covName)
  print(g_plot)
  return(g_plot)
}


##
# diagonstic tools for validating the model
#
##
diagonstic.geomm <- function(geomm, new.data = NULL, n.groups = 10){


  beta <- geomm$beta
  if(is.null(new.data)){
    y <- geomm$y
    X <- geomm$X
    offset <- geomm$offset
  }else{
    mf <- model.frame(formula = geomm$formula, data = new.data)
    y  <- as.matrix(model.extract(mf, "response"))
    X  <- as.matrix(model.matrix(attr(mf, "terms"), data = mf))
    offset <- model.offset(mf)
  }
  y_max <- max(y)
  u_y   <- unique(y)
  n_unique <- length(u_y)
  if(n_unique < 15){

    groups = sort(
      unique(
        round(c(0,quantile(unique(y), seq(0,1,length=n.groups))))))
  }else{

    groups <- sort(
      unique(
        round(c(0,quantile(y, seq(0,1,length=n.groups))))))

  }
  if(is.null(offset)){
    theta <- 1 + exp(X %*% beta)
  }else{
    theta <- 1 + exp(X %*% beta + offset)
  }

  n_dig <- max(nchar(y))
  string_y <- paste('%', n_dig,'d',sep="")
  cat(strrep('*',50 + 2*n_dig), '\n')
  cat( stringr::str_pad('Probability', 11 + 4 + 2*n_dig, 'right'), 'estimated', '(numerical)\n')
  for(i in 1:(length(groups)-1)){

    if(groups[i+1] - groups[i] == 1){
      index <- y==groups[i+1]
      P = dgeom(groups[i+1]-1, prob = 1/theta)
      string_1 <- paste('P(Y  =   ',sprintf(string_y,as.integer(groups[i+1])),' )  = ',sep="")
      string_1 <-  stringr::str_pad(string_1, nchar(string_1) + n_dig + 3, "right")
      mean_index <- mean(index)
      mean_P     <- mean(P)
      if(mean_index < 0.01){
        string_mean  <-sprintf('%0.2e',mean_index)
        string_P     <-sprintf('%0.2e',mean_P)
      }else{
        string_mean  <-sprintf('%2.4f',mean_index)
        string_P     <-sprintf('%0.4f',mean_P)
      }
      cat(string_1, string_P,' (',string_mean,')',' n = ',sum(index),'\n')
    }else{
      P_o = pgeom(groups[i+1]-1, prob = 1/theta)
      P_u = pgeom(groups[i]-1, prob = 1/theta)
      index = (groups[i] <y) * (y<=groups[i+1]) ==1

      mean_index <- mean(index)
      mean_P     <- mean(P_o-P_u)
      if(mean_index < 0.01){
        string_mean  <-sprintf('%0.2e',mean_index)
        string_P     <-sprintf('%0.2e',mean_P)
      }else{
        string_mean  <-sprintf('%2.4f',mean_index)
        string_P     <-sprintf('%0.4f',mean_P)
      }

      cat('P(Y in (',
          sprintf(string_y, as.integer(groups[i])),
          ',',
          sprintf(string_y, as.integer(groups[i+1])),
          ']) = ',string_P,' (',
          string_mean,')',' n = ',sum(index),'\n')
    }

  }
}




diagonstic.geomcm <- function(geomcm, new.data = NULL, K = NULL, n.groups = 10){


  beta <- geomcm$beta
  if(is.null(new.data)){
    y <- geomcm$y
    X <- geomcm$X
    offset <- geomcm$offset
    K <- geomcm$K
  }else{
    mf <- model.frame(formula = geomcm$formula, data = new.data)
    y  <- as.matrix(model.extract(mf, "response"))
    X  <- as.matrix(model.matrix(attr(mf, "terms"), data = mf))
    offset <- model.offset(mf)
  }
  y_max <- max(y)
  u_y   <- unique(y)
  n_unique <- length(u_y)
  if(n_unique < 15){

    groups = sort(
      unique(
        round(c(0,quantile(unique(y), seq(0,1,length=n.groups))))))
  }else{

    groups <- sort(
      unique(
        round(c(0,quantile(y, seq(0,1,length=n.groups))))))

  }
  theta <- rep(1, length(K))
  if(is.null(offset)){

    theta[K>2] <- 1 + exp(X %*% beta)
  }else{
    theta[K>2] <- 1 + exp(X %*% beta + offset)
  }

  n_dig <- max(nchar(y))
  string_y <- paste('%', n_dig,'d',sep="")
  cat(strrep('*',50 + 2*n_dig), '\n')
  cat(str_pad('Probability', 11 + 4 + 2*n_dig, 'right'), 'estimated', '(numerical)\n')

  P_K <- pgeom(K, prob = 1/theta)
  for(i in 1:(length(groups)-1)){

    if(groups[i+1] - groups[i] == 1){
      index <- y==groups[i+1]
      P <- rep(0, length(K))
      index_K = K >= groups[i+1]
      P[index_K] <- dgeom(groups[i+1]-1, prob = 1/theta[index_K])/P_K[index_K]
      string_1 <- paste('P(Y  =   ',sprintf(string_y,as.integer(groups[i+1])),' )  = ',sep="")
      string_1 <- stringr::str_pad(string_1, nchar(string_1) + n_dig + 3, "right")
      mean_index <- mean(index)
      mean_P     <- mean(P)
      if(mean_index < 0.01){
        string_mean  <-sprintf('%0.2e',mean_index)
        string_P     <-sprintf('%0.2e',mean_P)
      }else{
        string_mean  <-sprintf('%2.4f',mean_index)
        string_P     <-sprintf('%0.4f',mean_P)
      }
      cat(string_1, string_P,' (',string_mean,')',' n = ',sum(index),'\n')
    }else{
      P <- rep(0, length(K))
      index_K = K >= groups[i+1]
      P[index_K] = pgeom(groups[i+1]-1, prob = 1/theta[index_K])
      P[index_K] = P[index_K] - pgeom(groups[i]-1, prob = 1/theta[index_K])
      P[index_K] = P[index_K] / P_K[index_K]
      index = (groups[i] <y) * (y<=groups[i+1]) ==1

      mean_index <- mean(index)
      mean_P     <- mean(P)
      if(mean_index < 0.01){
        string_mean  <-sprintf('%0.2e',mean_index)
        string_P     <-sprintf('%0.2e',mean_P)
      }else{
        string_mean  <-sprintf('%2.4f',mean_index)
        string_P     <-sprintf('%0.4f',mean_P)
      }

      cat('P(Y in (',
          sprintf(string_y, as.integer(groups[i])),
          ',',
          sprintf(string_y, as.integer(groups[i+1])),
          ']) = ',string_P,' (',
          string_mean,')',' n = ',sum(index),'\n')
    }

  }
}
##
# prediction
##
predict.geomcm <- function(geomcmObj, newdata = NULL, response = 'link'){
  if(is.null(newdata)){
    X <- geomcmObj$X
    offset <- geomcmObj$offset
  }else{
    if(as.character(getResponseFormula(as.formula(geomcmObj$formula))[[2]] )%in% colnames(newdata) == FALSE)
      newdata$y <- rep(2,nrow(newdata))
    init_list <- init_geomm(geomcmObj$formula, newdata)
    X = init_list$X
    offset = init_list$offset
    if(is.null(offset))
      offset <- rep(0, length(newdata$y))
  }
  pA <- NULL
  if(response == 'link')
    pA <- (1+exp(X[,dimnames(geomcmObj$beta)[[1]]]%*%geomcmObj$beta + offset))

  return(pA)
}
predict.geomm <- function(geommObj, newdata = NULL, response = 'link'){

  if(is.null(newdata)){
    X <- geommObj$X
    offset <- geommObj$offset
  }else{
    if(as.character(getResponseFormula(as.formula(geommObj$formula))[[2]] )%in% colnames(newdata) == FALSE)
      newdata$y <- rep(2,nrow(newdata))
    init_list <- init_geomm(geommObj$formula, newdata)
    X = init_list$X
    offset = init_list$offset
    if(is.null(offset))
      offset <- rep(0, length(newdata$y))
  }
  pA <- NULL
  if(response == 'link')
    pA <-(1+exp(X%*%geommObj$beta + offset))

  return(pA)
}



