logit <-
  function (p) 
    log(p/(1 - p))

cox_first_degree = function(y, prob){
  # y = c(rep(0, 9), rep(1,1), rep(0,8), rep(1, 2), rep(0,7),rep(1,3), rep(0,6),rep(1,4), rep(0,5), rep(1,5), rep(0,4),rep(1,6),
  #       rep(0,3),rep(1,7), rep(0,2), rep(1,8), rep(0,1), rep(1,9))
  # # prob = c(rep(0.1, 10), rep(0.2, 10),rep(0.3, 10),rep(0.4, 10),rep(0.5, 10),rep(0.6, 10),rep(0.7, 10),rep(0.8, 10),rep(0.9, 10))
  # prob = c(rep(0.3, 10), rep(0.3, 10),rep(0.3, 10),rep(0.4, 10),rep(0.5, 10),rep(0.6, 10),rep(0.7, 10),rep(0.8, 10),rep(0.9, 10))
  # 
  dat <- data.frame(e = prob, o = y)
  dat$e[dat$e == 0] = 0.0000000001
  dat$e[dat$e == 1] = 0.9999999999
  dat$logite <- logit(dat$e)
  
  mfit = glm(formula = o~I(logite), 
             family = binomial(link = "logit"), dat)
  # browser()
  slope = mfit$coefficients[2]
  intercept = mfit$coefficients[1]
  return(list(slope = slope, intercept = intercept))
}


# xweight <- seq(-3, 3, 0.001)
# yweight <- predict(mfit, list(logite = xweight),type="response")
# quartz('test')
# plot(dat$e, dat$o, pch = 16, xlab = "WEIGHT (g)", ylab = "VS")
# lines(xweight, yweight)