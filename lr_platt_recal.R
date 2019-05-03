lr_platt_recal_func = function(model_fit, validate_data, test_res, test){
  # predict validation datasets's estimates with logistic regression
  val_estimates_norm = predict(model_fit, validate_data[,-length(validate_data)], type='response')
  train_re_mtx = cbind(y=validate_data$y, yhat=val_estimates_norm)
  
  # create calibration model
  calib.model <- glm(y~yhat, as.data.frame(train_re_mtx), family=binomial)
  ygrid_norm = as.data.frame(test_res)
  colnames(test_res) <- c("yhat")
  
  # recalibrate and measure on test set
  ygrid_cal = predict(calib.model, test_res, type='response')
 
  return(ygrid_cal)
}