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

cam_obj <- growthObjInit(reg_dat, "sales", "sales_prev", exitIndex=  "exit_term")
dist_list <- list(Above =     " log(sales/sales_prev) ~ 1  + ln_roa +  age + emp_prev  ",
                  Below =     " y ~ 1  + ln_roa +  age +  emp_prev ",
                  ProbEqual = " y ~ 1  + ln_roa + emp_prev ",
                  ProbAbove = " y ~ 1  + ln_roa + age +sales_prev ",
                  ProbExit  = " y ~ 1  + ln_roa + emp_prev + sales_prev ")
cam_obj <- growthEstimate(cam_obj,
                          formula = dist_list,
                          dist = list(Above = "gamma", Below = "betamc"),
                          difference=F)

logNorm_obj <- lm("log(sales/sales_prev) ~ 1 + ln_roa + age + emp_prev",
                  data = reg_dat)
cat('Above:\n')
print(summary(cam_obj$estimate$Above$model.fit))
cat('Below:\n')
print(summary(cam_obj$estimate$Below$model.fit))
cat('Prob above:\n')
print(summary(cam_obj$estimate$ProbAbove))
cat('log normal:\n')
print(summary(logNorm_obj))
