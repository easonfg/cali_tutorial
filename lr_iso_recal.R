
lr_iso_recal_func = function(model_fit, validate_data, test_res, test){
  # predict validation datasets's estimates with logistic regression
  val_estimates_norm = predict(model_fit, validate_data[,-length(validate_data)],  type='response')
  train_re_mtx = cbind(y=validate_data$y, yhat=val_estimates_norm)
  iso_train_mtx = train_re_mtx[order(train_re_mtx[,2]),]
  
  # create calibration model
  calib.model <- isoreg(iso_train_mtx[,2], iso_train_mtx[,1])
  stepf_data = cbind(calib.model$x, calib.model$yf)
  step_func = stepfun(stepf_data[,1], c(0,stepf_data[,2]))
  
  # recalibrate and measure on test set
  exp2_iso_recal <- step_func(test_res[,1])

  return(exp2_iso_recal)
}
