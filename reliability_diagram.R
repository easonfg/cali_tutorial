reliability_diagram = function(obs_rep, data_ls, stat_type,
                          title_ls, 
                          color_ls, title) {
  
  data_ls_len = length(data_ls)
  
  ### create bin averages
  for (data_i in 1:data_ls_len){
    temp_res = reliability_datapts(obs_rep, data_ls[[data_i]], bins = 10, stat_type = stat_type)
    assign(paste("recal_bins", data_i, sep=""),temp_res)
  }
  
  for (data_i in 1:data_ls_len){
    temp_res = melt(get(paste('recal_bins', data_i, sep='')), id="V2")
    temp_res[, "variable"] <- paste("Vol.x", data_i, sep='')
    assign(paste('melt', data_i, sep=''), temp_res)
  }
  
  data = melt1
  if (data_ls_len > 1){
    for (data_i in 2:data_ls_len){
      data = rbind(data, get(paste('melt', data_i, sep='')))
    }
  }

  line_plot = ggplot(data, aes(x=V2,  y=value, color=variable))   +  geom_point()+ geom_line() +
    scale_color_manual(labels = title_ls,
                       values = color_ls) +
    guides(color=guide_legend(" ")) + 
    ggtitle(paste(stat_type, 'Statistics')) +
    xlab("Mean Predicted") + ylab("Mean Observed") + 
    xlim(0, 1) + ylim(-0.05, 1.05)
  
  line_plot = line_plot + geom_abline(intercept = 0, slope = 1, color="black",
                                      linetype="dashed", size=1) 
  #remove background
  line_plot = line_plot + theme_bw()
  line_plot = line_plot + theme(legend.position="bottom")
  
  ## add data points
  for (data_i in 1:data_ls_len){
    temp_obs = obs_rep
    temp_obs[temp_obs==0] = -0.005 * data_i
    temp_obs[temp_obs==1] = 1 + 0.005 * data_i
    assign(paste('obs_rep_offset', data_i, sep=''), data.frame(cbind(data_ls[[data_i]], temp_obs)))
  }
  
  for (data_i in 1:data_ls_len){
    temp_res = melt(get(paste('obs_rep_offset', data_i, sep='')), id="temp_obs")
    temp_res[, "variable"] <- paste("Vol.x", data_i, sep='')
    assign(paste('obs_rep_offset', data_i, sep=''), temp_res)
  }
  
  data_points = obs_rep_offset1
  if (data_ls_len > 1){
    for (data_i in 2:data_ls_len){
      data_points = rbind(data_points, get(paste('obs_rep_offset', data_i, sep='')))
    }
  }
  
  line_plot = line_plot + geom_point(data = data_points, aes(x=data_points$value,  y=data_points$temp_obs, color=variable)) + geom_point(alpha=0.2)
  
  
  ggsave(title, device = 'png', width = 6, height = 6)
}

reliability_datapts <- function(obs, pred, bins=10, stat_type ='H') {
  min.pred <- min(pred)
  max.pred <- max(pred)
  min.max.diff <- max.pred - min.pred
  
  if (stat_type == 'H'){
    mtx = cbind(obs, pred)
    mtx = as.data.frame(mtx)
    mtx = mtx[order(mtx$pred),]
    res = data.frame(V1= numeric(0), V2 = numeric(0))
    split_mtx = split(mtx, cut(mtx$pred, seq(0,1,1/10), include.lowest=TRUE))

    for (i in 1:length(split_mtx)){
      col_mean = colMeans(split_mtx[[i]])
      if (sum(is.na(col_mean)) > 0) {
        next
      }
      res[i,] = col_mean
    }
    
  }else{
    ## C statistics, same number of instances in each bin
    mtx = cbind(obs, pred)
    mtx = as.data.frame(mtx)
    mtx = mtx[order(mtx$pred),]
    n <- length(pred)/10
    nr <- nrow(mtx)
    split_mtx = split(mtx, rep(1:ceiling(nr/n), each=n, length.out=nr))
    res = data.frame(V1= numeric(0), V2 = numeric(0))
    for (i in 1:length(split_mtx)){
      res[i,] = colMeans(split_mtx[[i]])
    }
  }
  
  return(res)
}