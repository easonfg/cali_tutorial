

###### plotting fitted lines ###### ###### ###### ###### ###### 
platt_func_draw = function(validate, val_estimates_norm){
  quartz(title="logistic") # creates a quartz window with title
  
  dat = as.data.frame(cbind(y=validate$y, yhat=val_estimates_norm[,1]))
  plot(dat$yhat,dat$y,
       xlim=c(-0.1, .4), xlab="Uncalibrated Estimates",ylab="Calibrated Estimates") # plot with body size on x-axis and survival (0 or 1) on y-axis
  logistic_model=glm(y~.,family=binomial,dat) # run a logistic regression model (in this case, generalized linear model with logit link). see ?glm
  
  curve(predict(logistic_model,data.frame(yhat=x),type="resp"),add=TRUE, col = "blue", lwd = 2) # draws a curve based on prediction from logistic regression model
  legend("left", c('Platt scaling fit','Data'),lty=c(1,NA),pch=c('','o'), col = c('blue', 'black'),bg='white')
  quartz.save('logfit.jpg', type = "jpg")
  
  
}