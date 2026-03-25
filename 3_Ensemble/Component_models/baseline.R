### ------------------------------------------------------------------------ ### 
### Baseline model 
### ------------------------------------------------------------------------ ### 

# Function taken from https://github.com/sbfnk/covid19.forecasts.uk 
null_model_forecast_quantiles <- function(values, horizon, truncation = 0)
{
  if (truncation > 0) {
    values <- head(values, -truncation)
  }
  
  tnorm_unc_fit <- function(x, mean, true) {
    return(-sum(dtnorm(x = true, mean = mean, sd = x, lower = 0, log = TRUE)))
  }
  
  quantiles <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
  ts <- tail(values, min(horizon, length(values)))
  
  if (length(ts) > 1) {
    ## null model
    null_ts <- rep(last(ts), horizon + truncation)
    
    interval <- c(0, max(abs(diff(ts))))
    
    if (diff(interval) > 0 ) {
      tnorm_sigma <-
        optimise(tnorm_unc_fit, interval = interval,
                 mean = head(ts, -1),
                 true = tail(ts, -1))
      
      quant <- round(qtnorm(p = quantiles, mean = last(ts),
                            sd = tnorm_sigma$minimum, lower = 0))
    } else {
      quant <- rep(last(ts), length(quantiles))
    }
  } else {
    quant <- rep(last(ts), length(quantiles))
  }
  
  names(quant) <- quantiles
  
  return(quant)
}


res <- tibble()

for (iRegion in seq_along(regions)) {
  
  myreg<-regions[iRegion]
  print(myreg)
  
  for (iDate in seq_along(dates_analysis)) {
    
    print(dates_analysis[iDate])
    
    dat_train<-data_rt %>%
      filter(date_analysis==dates_analysis[iDate],
             region==myreg)  %>%
      as.data.frame()
    
    # FIT MODEL
    quant <-
      null_model_forecast_quantiles(dat_train$iHosp_smooth, horizon = 14 +remove_last_n,
                                    truncation = 0)
    
    # FORMAT DATA
    names(quant)<-c(paste0("lower_",c(2,5,seq(10,90,by=10))),"point",paste0("upper_",rev(c(2,5,seq(10,90,by=10)))))
    res_tmp <- data.frame(
      var = "iHosp",
      region = myreg, 
      date = dates_analysis[iDate],
      prediction_horizon = seq(1,max_prediction_horizon)-remove_last_n,
      t(quant)
      
    )
    
    res <- rbind(res, res_tmp)
    
    
  }
  
}


res <- res %>% mutate(model = model)

write.csv2(res,file=paste0(path_results,model,"_period",period,".csv"),row.names=F)
