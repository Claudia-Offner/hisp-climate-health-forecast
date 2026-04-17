### ------------------------------------------------------------------------ ### 
### ARIMA1: Autoregressive integrated moving average model 
### ------------------------------------------------------------------------ ### 

require(forecast)

res<-tibble()

for (iLocation in seq_along(locations)) {
  
  myloc<-locations[iLocation]

  for (iDate in seq_along(dates_analysis)){
    
    print(paste(model, myloc, dates_analysis[iDate]))
    
    dat_train<-data_rt %>%
      filter(date_analysis==dates_analysis[iDate],
             location==myloc) 
    
    # FIT MODEL (w/ week and season effects)
    # mod<-auto.arima(ts(dat_train$outcome, frequency=7))
    dat_msts <- msts(dat_train$outcome, seasonal.periods = c(7, 365.25))
    mod <- tbats(dat_msts)
    proj<-forecast(mod,max_prediction_horizon)
    
    # FORMAT DATA
    out <- as.data.frame(forecast(mod,max_prediction_horizon,level=c(10,20,30,40,50,60,70,80,90,95,98)))
    out<-sapply(out,FUN=function(x){pmax(0,round(x))})
    alpha<-rev(c(2,5,10,20,30,40,50,60,70,80,90))
    colnames(out)<-c("point",paste0(c("lower_","upper_"),rep(alpha,each=2)))
    
    res_tmp<- data.frame(var="outcome",
                         location=myloc,
                         date=dates_analysis[iDate],
                         prediction_horizon=seq(1,max_prediction_horizon)-remove_last_n,
                         out)
    
    res<-rbind(res,res_tmp)
    
  }
}


res <- res %>% mutate(model=model)

write.csv(res,file=paste0(path_results,model,"_period",period,".csv"),row.names=F)
