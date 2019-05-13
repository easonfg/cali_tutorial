source('./platt_func_draw.R')

svm_platt_recal_func = function(model_fit, validate_data, test_res, test){
  # predict validation datasets's estimates with svm
  val_estimates = predict(model_fit, validate_data[,-length(validate_data)], probability=TRUE)
  val_estimates_norm = as.data.frame(attr(val_estimates, 'probabilities')[,2])
  train_re_mtx = cbind(y=as.numeric(validate_data$y), yhat=val_estimates_norm[,1])


  
  # create calibration model
  calib.model <- glm(y~yhat, as.data.frame(train_re_mtx), family=binomial)
  ygrid_norm = as.data.frame(test_res)
  colnames(ygrid_norm) <- c("yhat")
  
  # recalibrate and measure on test set
  ygrid_cal = predict(calib.model, ygrid_norm, type='response')
  
  ### platt fit func draw
  platt_func_draw(validate_data, val_estimates_norm)
  ###
  
  return(ygrid_cal)
}
