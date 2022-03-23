
logit <- function(x){1/(1+exp(-x))}

#'
#' simulate company data that mimics the article
#'
#'
#'
#'
generate.data <- function(n, size="small"){
  if(size=="small"){
    emp_prev <- rgeom(n, 1/2.06)+1

    emp_prev <- rgeom(n, 1/2.06)
    while(sum(emp_prev>8)>0){
      index =  emp_prev>8
      emp_prev[index] = rgeom(sum(index),1/2.06)
    }
    emp_prev <- emp_prev+1
    #regression on the other variables
    #model.age = glm.nb(age-2 ~ 1 + emp_prev, data=reg_dat)
    age = rnbinom(n, mu = exp(0.9397 ), size =  1.4462)  +2
    #lm_roa.model <- lm(ln_roa ~ emp_prev + age, dat=  reg_dat)
    ln_roa <-  0.1088506  + -0.0021723  * emp_prev + -0.0032622 * age + rnorm(n, sd = 0.1037)
    #assp.model <-  lm(assp_hist ~ emp_prev + age + ln_roa , dat=  reg_dat)
    assp_hist <- -0.24157 + 0.52207 * emp_prev + -0.20382 * age + 7.9261 * ln_roa + rnorm(n, sd = 6.837)

    reg_dat.simulate <- data.frame( emp_prev        = emp_prev,
                                    age             = age,
                                    ln_roa          = ln_roa,
                                    assp_hist       = assp_hist,
                                    assp_hist_over  = assp_hist,
                                    assp_hist_below = assp_hist)
    reg_dat.simulate$assp_hist_over[reg_dat.simulate$assp_hist_over<0] = 0
    reg_dat.simulate$assp_hist_below[reg_dat.simulate$assp_hist_below>0] = 0

    X = cbind(rep(1,n),
              reg_dat.simulate$ln_roa,
              reg_dat.simulate$age,
              reg_dat.simulate$assp_hist_over,
              reg_dat.simulate$assp_hist_below)
    ##
    # simulate model
    #
    ##

    #exit_term
    beta.exit <- c(2.14, 0.1,0, -0.04, 0.125)
    Emp_effect <- -0.5*emp_prev>1
    Exit <- runif(n) > logit(X%*%beta.exit +  Emp_effect)
    #equal
    beta.equal <- c(-1.41, -2.39, 0, 0.03, -0.03)
    Emp_effect <- -3 + 3*exp(1-emp_prev)
    Equal <- runif(n) > logit(X%*%beta.equal +  Emp_effect)
    #prob above
    beta.above <- c(0.23, 3.42, 0, -0.19, -0.22)
    Above <- runif(n) > logit(X%*%beta.above )
    # Above
    beta.Above <- c(-0.86, -1.71, -0.142, 0.06, -0.07)
    Above.val = rgeom(n, prob = 1/(1+exp(X%*%beta.Above+ log(emp_prev))))  + emp_prev
    # Below
    beta.Below <- c(-2.71, -4.40, -0.01, 0.33, -0.18)
    Belov.val = rgeomc(X, beta.Below, emp_prev, offset =log(emp_prev))

    reg_dat.simulate$exit <- Exit
    reg_dat.simulate$emp <- emp_prev*Equal +
                            (1-Equal) * Above * Above.val+
                            (1-Equal) * (1-Above) * Belov.val

  }else{

    emp_prev <- rgeom(n, 1/8.9)
    while(sum(emp_prev>40)>0){
      index =  emp_prev>40
      emp_prev[index] = rgeom(sum(index),1/8.9)
    }
    emp_prev <- emp_prev+10
    #model.age = glm.nb(age-2 ~ 1 + emp_prev, data=reg_dat)
    age = rnbinom(n, mu = exp(1.458010 + 0.003768*emp_prev  ), size =  1.821)  + 2
    #lm_roa.model <- lm(ln_roa ~ emp_prev + age, dat=  reg_dat)
    ln_roa <-  8.859e-02  + -8.679e-04  * emp_prev + 8.498e-05 * age + rnorm(n, sd = 0.09403)
    #assp.model <-  lm(assp_hist ~ emp_prev + age + ln_roa , dat=  reg_dat)
    assp_hist <- 2.19752 +-0.01359   * emp_prev + -0.09728  * age + 3.36240 * ln_roa + rnorm(n, sd = 6.498 )

    reg_dat.simulate <- data.frame( emp_prev        = emp_prev,
                                    age             = age,
                                    ln_roa          = ln_roa,
                                    assp_hist       = assp_hist,
                                    assp_hist_over  = assp_hist,
                                    assp_hist_below = assp_hist)
    reg_dat.simulate$assp_hist_over[reg_dat.simulate$assp_hist_over<0] = 0
    reg_dat.simulate$assp_hist_below[reg_dat.simulate$assp_hist_below>0] = 0

    X = cbind(rep(1,n),
              reg_dat.simulate$ln_roa,
              reg_dat.simulate$age,
              reg_dat.simulate$assp_hist_over,
              reg_dat.simulate$assp_hist_below)

    beta.exit <- c(-2.27, -3.74,0, 0.05, 0.015)
    Exit <- runif(n) > logit(X%*%beta.exit)
    #equal
    beta.equal <- c(-1.54, 0.44, 0, -0.01, 0.02)
    Equal <- runif(n) > logit(X%*%beta.equal)
    #prob above
    beta.above <- c(0.27, 3.52, 0, -0.09, -0.08)
    Above <- runif(n) > logit(X%*%beta.above )
    # Above
    beta.Above <- c(-1.65, -0.2, -0.03, -0.01, -0.02)
    Above.val = rgeom(n, prob = 1/(1+exp(X%*%beta.Above + log(emp_prev))))  + emp_prev
    # Below
    beta.Below <- c(-2.31, -0.69, -0.02, 0.11, -0.06)
    Belov.val =  rgeomc(X, beta.Below, emp_prev, offset = log(emp_prev))

    reg_dat.simulate$exit <- Exit
    reg_dat.simulate$emp <- emp_prev*Equal +
      (1-Equal) * Above * Above.val+
      (1-Equal) * (1-Above) * Belov.val

  }
  return(reg_dat.simulate)
}
