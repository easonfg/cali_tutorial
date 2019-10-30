ici = function(Y, P){
  # Y = c(rep(0, 9), rep(1,1), rep(0,8), rep(1, 2), rep(0,7),rep(1,3), rep(0,6),rep(1,4), rep(0,5), rep(1,5), rep(0,4),rep(1,6),
  #       rep(0,3),rep(1,7), rep(0,2), rep(1,8), rep(0,1), rep(1,9))
  # # P = c(rep(0.1, 10), rep(0.2, 10),rep(0.3, 10),rep(0.4, 10),rep(0.5, 10),rep(0.6, 10),rep(0.7, 10),rep(0.8, 10),rep(0.9, 10))
  # P = c(rep(0.3, 10), rep(0.3, 10),rep(0.3, 10),rep(0.4, 10),rep(0.5, 10),rep(0.6, 10),rep(0.7, 10),rep(0.8, 10),rep(0.9, 10))
  loess.calibrate <- loess(Y ~ P)
  # 
  # Estimate loess‐based smoothed calibration curve
  P.calibrate <- predict (loess.calibrate, newdata = P)
  
  # This is the point on the loess calibration curve corresponding to a given predicted probability.
  ICI <- mean (abs(P.calibrate - P))
  # browser()
  return(ICI)
  quartz('test')
  plot(P, P.calibrate)
  # plot(main="Loess Smoothing and Prediction", xlab="Date", ylab="Unemployment (Median)")
  lines(P.calibrate, x=P)
}