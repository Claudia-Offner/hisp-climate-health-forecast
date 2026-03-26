# Compute the ensemble forecasts for hospital admissions and project ICU and beds  


# Set path and parameters ---- 

# size of the window (in days) on which the paramaters for beds are estimated
length_calib<-5


# Load and combine projections from component models ----

proj_components<-tibble()
for (model in models){
  proj_tmp<-read.csv2(file=paste0("3_Ensemble/Results/",model,"_period",period,".csv"),stringsAsFactors =F)
  proj_components<-rbind.fill(proj_components,proj_tmp)
}

proj_components <- proj_components %>%
  mutate(date=as.Date(date)) 


# Compute ensemble average ---- 

ens<-proj_components %>%
  filter(model %in% models_ensemble) %>%
  group_by(var,region,date,prediction_horizon) %>%
  dplyr::summarise(point=mean(point,na.rm=T),
                   lower=mean(lower_5,na.rm=T),
                   upper=mean(upper_5,na.rm=T),
                   lower_2=mean(lower_2,na.rm=T),
                   upper_2=mean(upper_2,na.rm=T),
                   lower_5=mean(lower_5,na.rm=T),
                   upper_5=mean(upper_5,na.rm=T),
                   lower_10=mean(lower_10,na.rm=T),
                   upper_10=mean(upper_10,na.rm=T),
                   lower_20=mean(lower_20,na.rm=T),
                   upper_20=mean(upper_20,na.rm=T),
                   lower_30=mean(lower_30,na.rm=T),
                   upper_30=mean(upper_30,na.rm=T),
                   lower_40=mean(lower_40,na.rm=T),
                   upper_40=mean(upper_40,na.rm=T),
                   lower_50=mean(lower_50,na.rm=T),
                   upper_50=mean(upper_50,na.rm=T),
                   lower_60=mean(lower_60,na.rm=T),
                   upper_60=mean(upper_60,na.rm=T),
                   lower_70=mean(lower_70,na.rm=T),
                   upper_70=mean(upper_70,na.rm=T),
                   lower_80=mean(lower_80,na.rm=T),
                   upper_80=mean(upper_80,na.rm=T),
                   lower_90=mean(lower_90,na.rm=T),
                   upper_90=mean(upper_90,na.rm=T)) %>%
  mutate(model="Ensemble") 

proj<-bind_rows(proj_components,ens) %>%
  filter(date>=min(dates_analysis))  %>%
  filter(date<=max(dates_analysis)) %>%
  mutate(date_proj=date+prediction_horizon) 



# Project ICU and beds from iHosp Ensemble --------------------------------------------------------------------------------

simulbeds<-NULL


for (iRegion in seq_along(regions)){
  
  myregion<-regions[iRegion]
  
  
  for (iDate in seq_along(dates_analysis)) {
    
    # Point estimate ---------
    
    obs<- data_rt %>%
      filter(region==myregion) %>%
      dplyr::select(date,date_analysis,iHosp_smooth,iICU_smooth,inHosp_smooth,inICU_smooth,iHosp_smooth_lwr,iHosp_smooth_upr)%>%
      dplyr::rename("iHosp"="iHosp_smooth","iICU"="iICU_smooth","inHosp"="inHosp_smooth","inICU"="inICU_smooth") %>%
      filter(date_analysis==dates_analysis[iDate]) %>%
      filter(date<=dates_analysis[iDate]-remove_last_n) %>%
      drop_na(iHosp)
    
    model_sims<-obs %>%
      dplyr::select(date,iHosp)
    
    ens_tmp<-ens %>%
      ungroup() %>%
      filter(region==myregion & var=="iHosp" & date==dates_analysis[iDate]) %>%
      mutate(date=date+prediction_horizon) %>%
      dplyr::select(date,point,lower,upper) %>%
      pivot_longer(cols=c(point,lower,upper),names_to="type",values_to="iHosp")
    
    model_sims_proj<-model_sims %>%
      bind_rows(ens_tmp[ens_tmp$type=="point",c("date","iHosp")])
    
    # Calibration 
    start_fit_at <- max(obs$date)-length_calib+1
    n_to_forecast <- max_prediction_horizon
    # If < 0 coefficient of variations will be estimated
    fixed_cv_hosp <- 0.9
    fixed_cv_ICU <- 0.8
    estimates <- calibrateV2(obs, model_sims, start_fit_at,
                             fixed_cv_hosp = fixed_cv_hosp,
                             fixed_cv_ICU = fixed_cv_ICU,
                             n_to_forecast = n_to_forecast,
                             project_pICU = "const") # const
    estimates
    
    # Projections
    proj_point <- project(model_sims_proj, estimates$pICUs, estimates$mu_tauHosp_R,
                          estimates$cv_tauHosp_R, estimates$mu_tauICU_R,
                          estimates$cv_tauICU_R) %>%
      mutate(type="point")
    
    projections_point <- proj_point %>%
      gather(-date, -type, key = "var", value = "value") %>%
      pivot_wider(id_cols=c(date,var),names_from = "type",values_from="value") %>%
      filter(date > (dates_analysis[iDate]-remove_last_n) ) %>%
      mutate(date_proj=date) %>%
      mutate(date=dates_analysis[iDate]) %>%
      mutate(prediction_horizon=as.numeric(date_proj-date)) 
    
    
    # Confidence interval (bootstrap) ------------------
    
    res.boot<-NULL
    
    date.first= NULL # beginning of smoothing - if null the smallest date in data
    date.last=NULL   # end of smoothng - if null the largest date in data
    #parameter for the smoother
    bwd=8 # initial bandwidth to compute daily pattern; the data should have at least 2*bwd+4*7 days
    weekday=T # correct for weekly pattern 
    d=3 # local polynomial of degree 3
    dZ=3 # make quadratic for the last segment (linear=2, constant=1)
    # grid.bw=NULL ## bandwidth fof final smoothing, if a list then CV is computed for each value and the best is used
    # grid.d.sig=NULL ## parameter for least revision, value 0 is often best for Cubic/Quadratic, but otherwise choose best CV from list
    grid.bw=15
    grid.d.sig=0
    area<-1
    name<-myregion
    nSim<-100
    
    obs<- data_rt %>%
      filter(region==myregion) %>%
      dplyr::select(date,date_analysis,iHosp,iICU,inHosp,inICU)%>%
      filter(date_analysis==dates_analysis[iDate]) %>%
      filter(date<=dates_analysis[iDate]-remove_last_n) %>%
      drop_na(iHosp) %>%
      mutate(inICU = ifelse(is.na(inICU), 0, inICU))
    
    dat_train<-obs %>%
      mutate(day=date, name=name,area=area) %>%
      as.data.frame() 
    
    # bootstrapped smoothed trajectories 
    for (myvar in c("iHosp","iICU","inHosp","inICU")) {
      
      log.incid<-paste0("log.",myvar)
      dat_train[,paste0("log.",myvar)]<-log(dat_train[,myvar]+1)
      
      first.date = min(dat_train$day)
      last.date = max(dat_train$day)
      n.weeks = floor(as.numeric(last.date-first.date)/7)
      
      smoothed.data = real.time.smoother(data=dat_train,log.incid= paste0("log.",myvar), n.weeks=n.weeks,
                                         weekday=weekday,bwd=bwd, d=d, dZ=dZ,
                                         grid.bw=grid.bw, grid.d.sig=grid.d.sig)
      res.boot_tmp=boot.traj.LPR(smoothed.data,d=3,dZ=3,R=nSim)
      res.boot_tmp$y<-exp(res.boot_tmp$y)-1
      res.boot_tmp$var<-myvar
      res.boot<-bind_rows(res.boot,res.boot_tmp)
      
    }
    
    projections_boot<-NULL
    
    for (k in 1:nSim){
      
      obs<-res.boot %>%
        filter(i==k) %>% 
        dplyr::select(y,day,var) %>%
        pivot_wider(id_cols=day,names_from=var,values_from=y)%>%
        dplyr::rename(date=day) %>%
        drop_na()
      
      model_sims<-obs %>%
        dplyr::select(date,iHosp)
      
      ens_tmp<-ens %>%
        ungroup() %>%
        filter(region==myregion & var=="iHosp" & date==dates_analysis[iDate]) %>%
        mutate(date=date+prediction_horizon) %>%
        dplyr::select(date,point,lower,upper) %>%
        pivot_longer(cols=c(point,lower,upper),names_to="type",values_to="iHosp")
      
      model_sims_proj<-model_sims %>%
        bind_rows(ens_tmp[ens_tmp$type=="point",c("date","iHosp")])
      
      model_sims_proj_lwr<-model_sims %>% 
        bind_rows(ens_tmp[ens_tmp$type=="lower",c("date","iHosp")])
      
      model_sims_proj_upr<-model_sims %>% 
        bind_rows(ens_tmp[ens_tmp$type=="upper",c("date","iHosp")])
      
      # Calibration
      start_fit_at <- max(obs$date)-length_calib+1
      n_to_forecast <- max_prediction_horizon
      # If < 0 coefficient of variations will be estimated
      fixed_cv_hosp <- 0.9
      fixed_cv_ICU <- 0.8
      nrow(obs)
      nrow(model_sims)
      

      estimates <- calibrateV2(obs, model_sims, start_fit_at,
                               fixed_cv_hosp = fixed_cv_hosp,
                               fixed_cv_ICU = fixed_cv_ICU,
                               n_to_forecast = n_to_forecast,
                               project_pICU = "const") # const
      estimates
      
      # Projections
      proj_point <- project(model_sims_proj, estimates$pICUs, estimates$mu_tauHosp_R,
                            estimates$cv_tauHosp_R, estimates$mu_tauICU_R,
                            estimates$cv_tauICU_R) %>%
        mutate(type="point")
      proj_lwr <- project(model_sims_proj_lwr, estimates$pICUs, estimates$mu_tauHosp_R,
                          estimates$cv_tauHosp_R, estimates$mu_tauICU_R,
                          estimates$cv_tauICU_R) %>%
        mutate(type="lower")
      proj_upr <- project(model_sims_proj_upr, estimates$pICUs, estimates$mu_tauHosp_R,
                          estimates$cv_tauHosp_R, estimates$mu_tauICU_R,
                          estimates$cv_tauICU_R) %>%
        mutate(type="upper")
      
      projections_tmp <- proj_point %>%
        bind_rows(proj_lwr) %>%
        bind_rows(proj_upr) %>%
        gather(-date, -type, key = "var", value = "value") %>%
        pivot_wider(id_cols=c(date,var),names_from = "type",values_from="value") %>%
        filter(date > (dates_analysis[iDate]-remove_last_n) ) %>%
        mutate(date_proj=date) %>%
        mutate(date=dates_analysis[iDate]) %>%
        mutate(prediction_horizon=as.numeric(date_proj-date)) %>%
        mutate(i=k)
      
      projections_boot<-bind_rows(projections_boot,projections_tmp)
      
    }
    
    projections<-projections_boot %>%
      group_by(date,date_proj,prediction_horizon,var) %>%
      dplyr::summarise(
        lower=min(lower),
        upper=max(upper)) %>%
      left_join(projections_point) %>%
      mutate(region=myregion,model="Ensemble")
    
    simulbeds<-bind_rows(simulbeds,projections)
    
  }
}



proj<-bind_rows(proj,simulbeds)

saveRDS(proj,file=paste0(path_results,"Ensemble.rds"))

