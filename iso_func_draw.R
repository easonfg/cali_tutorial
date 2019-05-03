iso_func_draw = function(validate, val_estimates_norm, calib.model){
  quartz(title="iso") # creates a quartz window with title
  
  dat = as.data.frame(cbind(y=validate$y, yhat=val_estimates_norm[,1]))
  
  iso_reg_pts = as.data.frame(cbind(calib.model$x, calib.model$yf))
  sf1 = stepfun(iso_reg_pts$V1, c(0,iso_reg_pts$V2))
  plot.new()
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
  
  plot(dat$yhat,dat$y,
       xlim=c(-0.1, .4), xlab="Uncalibrated Estimates",ylab="Calibrated Estimates")
  
  lines(sf1, do.points = FALSE, col = "blue", lwd = 2)
  legend("left", c('Isotonic fit','Data'),lty=c(1,NA),pch=c('','o'), col = c('blue', 'black'),bg='white')
  quartz.save('isofit.jpg', type = "jpg")
}