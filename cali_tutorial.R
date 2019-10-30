rm(list=ls())

library(pROC)
library(ResourceSelection)
# library(randomForest)
library(e1071)
library(generalhoslem)
library(ggplot2)
library(reshape2)
library(plyr)
library(Hmisc)
library(rms)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source('./svm_platt_recal.R')
source('./svm_iso_recal.R')
source('./lr_platt_recal.R')
source('./lr_iso_recal.R')
source('./hosmer_lemeshow.R')
source('./Spiegelhalter_z.R')
source('./reliability_diagram.R')
source('./mce_ece.R')
source('./cox.R')
source('./ici.R')

set.seed(612)


  #####################

#### normalize to 0 and 1
vert_norm <- function(dataframe){
  dataframe = as.data.frame(dataframe)
  df_min = apply(dataframe, 2, min)
  df_max = apply(dataframe, 2, max)
  
  #matrix minus min of each col
  temp_df = t(t(dataframe)-apply(dataframe,2,min))
  #divided by max - min of each row
  return_df = t(t(temp_df) / (df_max - df_min))
  return(as.data.frame(return_df))
}
######


pipe_run = function(ninstances){

      
  ####### create 23 coefficients  ####### 
  iter <- 20
  out <- matrix(NA, nrow=iter, ncol=ninstances)
  counter = 1
  
  
  for (i in seq(from=3, to=60, by=3)){
    res = rbinom(ninstances, 1, i/100)
    out[counter,] = c(res)
    counter = counter + 1
  }
  
  out = rbind(out, rnorm(ninstances, 0.5, 0.1))
  out = rbind(out, rnorm(ninstances, 0.5, 0.5))
  out = rbind(out, rnorm(ninstances, 0.5, 1))
  
  
  B = rep(c(1.0, -0.5), times = c(1,1), length.out = 20)
  B = c(B, 1.0, -1.0, 1.0)
  
  logit_y = t(out)%*%B - 4.60
  prob = exp(logit_y)/(1+exp(logit_y))
  ####### create 23 coefficients  ####### 
  
  ##### use uniform distributions to create labels ####
  uni_dist = runif(ninstances, min = 0, max = 1)
  labels = as.integer(uni_dist < prob)
  
  
  dev_prob = rep(prob)
  dev_prob[dev_prob<0.027] = dev_prob[dev_prob<0.027]
  dev_prob[dev_prob>0.127] = dev_prob[dev_prob>0.127]

  dev_labels = as.integer(uni_dist < dev_prob)
  
  ##### use uniform distributions to create labels ####
  
  # combine coefficients and labels  
  perfect = as.data.frame(t(out))
  perfect = cbind(perfect, y=dev_labels)
  perfect = vert_norm(perfect)
  

  #### set training to be 50% of total number of samples and validation and test set to be 25%
  ### create train, validate, and test datasets
  smp_size <- floor(0.50 * nrow(perfect))
  val_smp_size = floor(0.25*nrow(perfect))

  train_ind <- sample(seq_len(nrow(perfect)), size = smp_size)
  train <- perfect[train_ind, ]
  
  val_test <- perfect[-train_ind, ]
  validate_ind <- sample(seq_len(nrow(val_test)), size = val_smp_size)
  validate <- val_test[validate_ind, ]
  test <- val_test[-validate_ind, ]
  #######
  

  
  ###
  
  ##### logistic regression ####
  ### train and predict  ### 
  #train with 'train' dataset and and predict 'test' dataset with logistic regressions 
  log_mod <- glm(y~., family=binomial(link='logit'),data=train)
  log_pred = as.data.frame(predict(log_mod, test[,-length(test)], type='response'))

  ### get auc  ### 
  pr <- ROCR::prediction(log_pred, test$y)
  auc_log <- ROCR::performance(pr, measure = "auc")
  auc_log <- auc_log@y.values[[1]]
  
  ### recalibration ###
  # lr recalibration with  platt
  lr_platt_recal = lr_platt_recal_func(log_mod, validate, log_pred, test)
  ## lr recalibration with  isotonic regression
  lr_iso_recal = lr_iso_recal_func(log_mod, validate, log_pred, test)
  ## recalibration ###
  

  ###### svm #######
  ### train and predict  ### 
  svmfit = svm(factor(y) ~ ., data = train, kernel = "linear", cost = 10, scale = FALSE, probability=TRUE)
  ygrid = predict(svmfit, test[,-length(test)], probability=TRUE)
  ygrid_norm = as.data.frame(attr(ygrid, 'probabilities')[,2])
  
  ### get auc  ### 
  pr <- ROCR::prediction(ygrid_norm, test$y)
  auc_svmnorm <- ROCR::performance(pr, measure = "auc")
  auc_svmnorm <- auc_svmnorm@y.values[[1]]
  
  ## svm recalibartion with  platt scaling
  svm_platt_recal = svm_platt_recal_func(svmfit, validate, ygrid_norm, test)
  ## svm recalibartion with  isotonic regression
  svm_iso_recal = svm_iso_recal_func(svmfit, validate, ygrid_norm, test)
  ###### svm #######
  
  ##average absolute difference ##
  print('log org avg absolute diff')
  print(mean(abs(test$y-log_pred[,1])))
  print('svm org avg aboslute diff')
  print(mean(abs(test$y-ygrid_norm[,1])))
  
  print('log platt avg absolute diff')
  print(mean(abs(test$y-lr_platt_recal)))
  print('svm platt avg absolute diff')
  print(mean(abs(test$y-svm_platt_recal)))
  
  print('log iso avg absolute diff')
  print(mean(abs(test$y-lr_iso_recal)))
  print('svm iso avg absolute diff')
  print(mean(abs(test$y-svm_iso_recal)))
  ##average absolute difference##
  
  
  ##calibration-in-the-large##
  print('log org cal in large')
  print(mean(test$y)/mean(log_pred[,1]))
  print('svm org cal in large')
  print(mean(test$y)/mean(ygrid_norm[,1]))
  
  print('log platt cal in large')
  print(mean(test$y)/mean(lr_platt_recal))
  print('svm platt cal in large')
  print(mean(test$y)/mean(svm_platt_recal))
  
  print('log iso cal in large')
  print(mean(test$y)/mean(lr_iso_recal))
  print('svm iso cal in large')
  print(mean(test$y)/mean(svm_iso_recal))
  ##calibration-in-the-large##
  
  ### measuring calibration ###
  ### hosmer lemeshow test ###
  print('hosmer lemeshow test: lr org C')
  hosmer_lemeshow(test$y, log_pred[,1], 10, 'C')
  print('hosmer lemeshow test: svm org C')
  hosmer_lemeshow(test$y, ygrid_norm[,1], 10, 'C')
  print('hosmer lemeshow test: lr platt C')
  hosmer_lemeshow(test$y, lr_platt_recal, 10, 'C')
  print('hosmer lemeshow test: lr iso C')
  hosmer_lemeshow(test$y, lr_iso_recal, 10, 'C')
  print('hosmer lemeshow test: svm platt C')
  hosmer_lemeshow(test$y, svm_platt_recal, 10, 'C')
  print('hosmer lemeshow test: svm iso C')
  hosmer_lemeshow(test$y, svm_iso_recal, 10, 'C')

  print('hosmer lemeshow test: lr org H')
  hosmer_lemeshow(test$y, log_pred[,1], 10, 'H')
  print('hosmer lemeshow test: svm org H')
  hosmer_lemeshow(test$y, ygrid_norm[,1], 10, 'H')
  print('hosmer lemeshow test: lr platt H')
  hosmer_lemeshow(test$y, lr_platt_recal, 10, 'H')
  print('hosmer lemeshow test: lr iso H')
  hosmer_lemeshow(test$y, lr_iso_recal, 10, 'H')
  print('hosmer lemeshow test: svm platt H')
  hosmer_lemeshow(test$y, svm_platt_recal, 10, 'H')
  print('hosmer lemeshow test: svm iso H')
  hosmer_lemeshow(test$y, svm_iso_recal, 10, 'H')

  ### hosmer lemeshow test ###
  
  ### Brier score ###
  cat('Brier score: svm org', val.prob(ygrid_norm[,1],test$y)['Brier'], '\n')
  cat('Brier score: lr org', val.prob(log_pred[,1],test$y)['Brier'], '\n')
  cat('Brier score: svm platt', val.prob(svm_platt_recal,test$y)['Brier'], '\n')
  cat('Brier score: lr platt', val.prob(lr_platt_recal,test$y)['Brier'], '\n')
  cat('Brier score: svm iso', val.prob(svm_iso_recal,test$y)['Brier'], '\n')
  cat('Brier score: lr iso', val.prob(lr_iso_recal,test$y)['Brier'], '\n')
  ### Brier score ###
  
  ### Spiegelhalter z test ###
  print('Spiegelhalter z test: lr org')
  Spiegelhalter_z(test$y, log_pred[,1])
  print('Spiegelhalter z test: svm org')
  Spiegelhalter_z(test$y, ygrid_norm[,1])
  print('Spiegelhalter z test: lr platt')
  Spiegelhalter_z(test$y, lr_platt_recal)
  print('Spiegelhalter z test: lr iso')
  Spiegelhalter_z(test$y, lr_iso_recal)
  print('Spiegelhalter z test: svm platt')
  Spiegelhalter_z(test$y, svm_platt_recal)
  print('Spiegelhalter z test: svm iso')
  Spiegelhalter_z(test$y, svm_iso_recal)

  ### Spiegelhalter z test ###
  
  #### mce ece ####
  print('lr org')
  print(ece_mce(test$y, log_pred[,1], 10, 'C'))
  
  print('svm org')
  print(ece_mce(test$y, ygrid_norm[,1], 10, 'C'))
  
  print('lr platt')
  print(ece_mce(test$y, lr_platt_recal, 10, 'C'))
  
  print('lr iso')
  print(ece_mce(test$y, lr_iso_recal, 10, 'C'))
  
  print('svm platt')
  print(ece_mce(test$y, svm_platt_recal, 10, 'C'))
  
  print('svm iso')
  print(ece_mce(test$y, svm_iso_recal, 10, 'C'))
  
  #### mce ece ####
  
  
  
  #### cox ####
  print('cox: lr org')
  print(cox_first_degree(test$y, log_pred[,1]))
  print('cox: svm org')
  print(cox_first_degree(test$y, ygrid_norm[,1]))
  print('cox: lr platt')
  print(cox_first_degree(test$y, lr_platt_recal))
  print('cox: lr iso')
  print(cox_first_degree(test$y, lr_iso_recal))
  print('cox: svm platt')
  print(cox_first_degree(test$y, svm_platt_recal))
  print('cox: svm iso')
  print(cox_first_degree(test$y, svm_iso_recal))
  #### cox ####
  
  #### ici ####
  print('ici: lr org')
  print(ici(test$y, log_pred[,1]))
  print('ici: svm org')
  print(ici(test$y, ygrid_norm[,1]))
  print('ici: lr platt')
  print(ici(test$y, lr_platt_recal))
  print('ici: lr iso')
  print(ici(test$y, lr_iso_recal))
  print('ici: svm platt')
  print(ici(test$y, svm_platt_recal))
  print('ici: svm iso')
  print(ici(test$y, svm_iso_recal))
  #### ici ####
  
  #### reliability diagrams  ####
  ###original LR and SVM with C statistics
  reliability_diagram(as.vector(test$y),
                 list(as.vector(log_pred[,1]),
                 as.vector(ygrid_norm[,1])),
                 'C',
                 c('LR Original', 'SVM Original'),
                 c('blue', 'red'),
                 'LR SVM original C')

  #original LR and SVM with H statistics
  reliability_diagram(as.vector(test$y),
                      list(as.vector(log_pred[,1]),
                           as.vector(ygrid_norm[,1])),
                      'H',
                      c('LR Original', 'SVM Original'),
                      c('blue', 'red'),
                      'LR SVM original H')

  #Platt and iso recalibrated LR with C statistics
  reliability_diagram(as.vector(test$y),
                list(as.vector(log_pred[,1]),
                as.vector(lr_platt_recal),
                as.vector(lr_iso_recal)),
                'C',
                c('LR Original', 'LR Platt scaling Recalibration', 'LR Isotonic regression Recalibration'),
                c('blue', 'green', 'orange'),
                'LR Platt Iso C')

  #Platt and iso recalibrated LR with H statistics
  reliability_diagram(as.vector(test$y),
                list(as.vector(log_pred[,1]),
                as.vector(lr_platt_recal),
                as.vector(lr_iso_recal)),
                'H',
                c('LR Original', 'LR Platt scaling Recalibration', 'LR Isotonic regression Recalibration'),
                c('blue', 'green', 'orange'),
                'LR Platt Iso H')

  #Platt and iso  recalibrated SVM with C statistics
  reliability_diagram(as.vector(test$y),
                 list(as.vector(ygrid_norm[,1]),
                 as.vector(svm_platt_recal),
                 as.vector(svm_iso_recal)),
                 'C',
                 c('SVM Original', 'SVM Platt scaling Recalibration', 'SVM Isotonic regression Recal'),
                 c('red', 'green', 'orange'),
                 'SVM Platt Iso C')

  #Platt and iso  recalibrated SVM with H statistics
  reliability_diagram(as.vector(test$y),
                 list(as.vector(ygrid_norm[,1]),
                 as.vector(svm_platt_recal),
                 as.vector(svm_iso_recal)),
                 'H',
                 c('SVM Original', 'SVM Platt scaling Recalibration', 'SVM Isotonic regression Recal'),
                 c('red', 'green', 'orange'),
                 'SVM Platt Iso H')
  ### measuring calibration ###
  
  
}



n5000res = pipe_run(ninstances = 5000)




