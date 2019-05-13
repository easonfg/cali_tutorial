hosmer_lemeshow = function(y, prob, g, stat_type){
  mtx = cbind(y, y_not = 1- y, prob, prob_not = 1-prob)
  mtx = as.data.frame(mtx)
  mtx = mtx[order(mtx$prob),]
  n <- length(prob)/g
  nr <- nrow(mtx)
  
  ## C statistics, same number of instances in each bin
  if (stat_type == 'C'){
    split_mtx = split(mtx, rep(1:ceiling(nr/n), each=n, length.out=nr))
  }else{ ### H statistics, equal intervals
    split_mtx = split(mtx, cut(mtx$prob, seq(0,1,1/g), include.lowest=TRUE))
    split_mtx = split_mtx[sapply(split_mtx, nrow)>0]
    ###
  }
  
  H_stat = 0
  for (i in 1:length(split_mtx)){
    obs = sum(split_mtx[[i]]$y == 1)
    exp = sum(split_mtx[[i]]$prob)
    obs_not = sum(split_mtx[[i]]$y == 0)
    exp_not = sum(split_mtx[[i]]$prob_not)
    
    if (exp == 0 || exp_not == 0){
      next
    }
    
    bin_sum = ((obs - exp)**2)/exp + ((obs_not - exp_not)**2)/exp_not

    H_stat = H_stat + bin_sum
  }
  PVAL = 1 - pchisq(H_stat, g - 2)
  
  cat('PVALUE', PVAL, '\n')
  cat('stat', H_stat, '\n')
  return(PVAL)
}