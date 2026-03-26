### ------------------------------------------------------------------------ ### 
### GAM2: Multivariate Generalized additive model
### ------------------------------------------------------------------------ ### 

library(mgcv)

if (period==1) {
  covar_all<-c('t2m_smooth', 'tp_smooth')
} else if (period==2){
  covar_all<-c('t2m_smooth', 'tp_smooth')
}

res<-tibble()

for (iDate in seq_along(dates_analysis)){
  
  print(dates_analysis[iDate]) 
  
  dat_train<-data_rt %>%
    filter(date_analysis==dates_analysis[iDate])  %>% 
    mutate(iHosp_r_covar=iHosp_r)
  
  # Estimate les best lags 
  res_lag<-tibble()
  list_col<-sort(covar_all)
  
  for (iCol in seq_along(list_col)){
    
    mycol<-list_col[iCol]
    
    for (lag in 1:12){
      
      sub1<-dat_train
      sub2<-sub1
      sub2$date<-sub2$date+lag
      
      sub<-left_join(sub1[,c("date","region","iHosp_r")],sub2[,c("date","region",mycol)], by = c("date", "region"))
      sub<-sub %>% drop_na(iHosp_r,mycol) %>%
        as.data.frame()
      cor<-cor.test(y=sub$iHosp_r,x=sub[,mycol],method="pearson")$estimate
      
      res_lag_tmp<-data.frame(pearson=cor,region="regions",lag=lag,var=mycol)
      res_lag<-bind_rows(res_lag,res_lag_tmp)
      
    }
  }
  
  
  best_lag<-res_lag %>% 
    group_by(var) %>%
    slice(which.max(abs(pearson)))
  if ("iHosp_r_covar" %in% list_col){
    best_lag$lag[best_lag$var=="iHosp_r_covar"]<-1
  }
  
  # Prepare lagged data
  list_lag<-best_lag$lag
  
  sub<-dat_train[,c("date","region","iHosp_r")]
  sub1<-dat_train[,c("date","region","iHosp_r",list_col)]
  
  for (i in 1:length(list_lag)){
    
    sub2<-sub1
    sub2$date<-sub2$date+list_lag[i]
    
    sub<-left_join(sub,sub2[,c("date","region",list_col[i])],by=c("region"="region","date"="date"))
    
  }
  
  sub<-sub %>% drop_na(iHosp_r,all_of(list_col)) %>%
    dplyr::select(-date,-region)
  
  # fit model on all regions 
  f<-as.formula(paste0("iHosp_r~",paste0("s(",list_col,")",collapse="+")))
  mod<-gam(f,data=sub)
  
  
  # predictions by region
  for (iRegion in seq_along(regions)) {
    
    myreg<-regions[iRegion]
    
    dat_train_reg<-dat_train %>%
      filter(region==myreg)
    
    # build new data to predict
    newdata<-data.frame(seq.Date(dates_analysis[iDate]-remove_last_n+1,dates_analysis[iDate]+14,by="day"),NA)
    names(newdata)<-c("date","iHosp_r")
    
    sub1<-dat_train_reg[,c("date","iHosp_r",list_col)]
    
    for (i in 1:length(list_lag)){
      
      sub2<-sub1
      sub2$date<-sub2$date+list_lag[i]
      
      newdata<-left_join(newdata,sub2[,c("date",list_col[i])],by=c("date"="date"))
      
      # fill NA by last available value
      newdata[is.na( newdata[,ncol(newdata)]),ncol(newdata)]<- last(newdata[!is.na( newdata[,ncol(newdata)]),ncol(newdata)])
      
    }
    
    # Predict growth rate
    proj<-predict(mod,newdata,se.fit=T)
    
    alpha<-c(2,5,10,20,30,40,50,60,70,80,90)/100
    out<-as.data.frame(proj$fit)
    
    for (iAlpha in seq_along(alpha)){ 
      
      lower<-proj$fit+qnorm(alpha[iAlpha]/2)*proj$se.fit
      upper<-proj$fit+qnorm(1-alpha[iAlpha]/2)*proj$se.fit
      
      out<-cbind(out,cbind(lower,upper))
      
    }
    
    
    # Predict hospitalizations 
    pred<-array(data=NA,dim=c(max_prediction_horizon,length(alpha)*2+1))
    pred[,1]<-dat_train_reg$iHosp_smooth[nrow(dat_train_reg)]*exp(cumsum(out[,1]))
    
    for (iAlpha in seq_along(alpha)){ 
      
      pred_tmp_lower<-dat_train_reg$iHosp_smooth_lwr[nrow(dat_train_reg)]*exp(cumsum(out[,iAlpha*2]))
      pred_tmp_upper<-dat_train_reg$iHosp_smooth_upr[nrow(dat_train_reg)]*exp(cumsum(out[,iAlpha*2+1]))
      pred[,c(2*iAlpha,2*iAlpha+1)]<-cbind(pred_tmp_lower,pred_tmp_upper)
      
    }
    
    pred<-sapply(as.data.frame(pred),FUN=function(x){pmax(0,round(x))})
    colnames(pred)<-c("point",paste0(c("lower_","upper_"),rep(alpha*100,each=2)))
    
    res_tmp<- data.frame(var="iHosp",
                         region=myreg,
                         date=dates_analysis[iDate],
                         prediction_horizon=seq(1,max_prediction_horizon)-remove_last_n,
                         pred)
    
    res<-rbind(res,res_tmp)
    
  }
}

res<-res%>%
  mutate(model=model)

write.csv2(res,file=paste0(path_results,model,"_period",period,".csv"),row.names=F)
