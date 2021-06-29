library(cams)
n <- 10000
size = "medium"
if(size== "small"){
  trunc_emp_prev=Inf
}else{
  name='medium'
  trunc_emp_prev = 15
}

data <- generate.data(n,size)
data$emp_prev_trunc <- data$emp_prev
data$emp_prev_trunc[data$emp_prev_trunc > trunc_emp_prev] = trunc_emp_prev
data$emp_prev_trunc <- factor(data$emp_prev_trunc)

cam_obj <- growthObjInit(data, "emp", "emp_prev", exitIndex=  "exit")
dist_list <- list(Above =     " y ~ 1  + ln_roa +  age + assp_hist_over + assp_hist_below +  offset(log(emp_prev))",
                  Below =     " y ~ 1  + ln_roa +  age + assp_hist_over + assp_hist_below +  offset(log(emp_prev))",
                  ProbEqual = " y ~ 1  + ln_roa + assp_hist_over + assp_hist_below + emp_prev_trunc ",
                  ProbAbove = " y ~ 1  + ln_roa + assp_hist_over + assp_hist_below ",
                  ProbExit  = " y ~ 1  + ln_roa + assp_hist_over + assp_hist_below + emp_prev_trunc ")
cam_obj <- growthEstimate(cam_obj, formula = dist_list)
logNorm_obj <- lm("log(emp) ~ -1 + assp_hist_over + assp_hist_below + emp_prev_trunc + offset(log(emp_prev))", data = data)
