### ------------------------------------------------------------------------ ### 
### GAM1: Generalized additive model
### ------------------------------------------------------------------------ ### 

res<-tibble()


for (iRegion in seq_along(regions)) {
  
  myreg<-regions[iRegion]
  print(myreg)
  
  for (iDate in seq_along(dates_analysis)){
    
    print(dates_analysis[iDate]) 
    
    dat_train<-data_rt %>%
      filter(date_analysis==dates_analysis[iDate],
             region==myreg) 
    dat_train$time<-1:nrow(dat_train)
    
    # FIT MODEL
    mod<-gam(iHosp_smooth~s(time,bs="ps",k=round(nrow(dat_train)/10)),data=dat_train)
    proj<-predict(mod,data.frame(time=seq(nrow(dat_train)+1,nrow(dat_train)+max_prediction_horizon)),se.fit=T)
    
    # FORMAT DATA
    alpha<-c(2,5,10,20,30,40,50,60,70,80,90)/100
    
    out<-as.data.frame(proj$fit)
    
    for (iAlpha in seq_along(alpha)){ 
      
      lower<-proj$fit+qnorm(alpha[iAlpha]/2)*proj$se.fit
      upper<-proj$fit+qnorm(1-alpha[iAlpha]/2)*proj$se.fit
      
      out<-cbind(out,cbind(lower,upper))
      
    }
    
    
    out<-sapply(out,FUN=function(x){pmax(0,round(x))})
    colnames(out)<-c("point",paste0(c("lower_","upper_"),rep(alpha*100,each=2)))
    
    
    res_tmp<- data.frame(var="iHosp",
                         region=myreg,
                         date=dates_analysis[iDate],
                         prediction_horizon=seq(1,max_prediction_horizon)-remove_last_n,
                         out)
    
    res<-rbind(res,res_tmp)
    
  }
  
}


res <- res %>% mutate(model=model)

write.csv2(res,file=paste0(path_results,model,"_period",period,".csv"),row.names=F)
